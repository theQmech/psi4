#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
Module of helper functions for ccresponse distributed property calculations.
Defines functions for interacting with the database created by the run_XXX
driver function.

Properties that are able to use this module should be added to
the registered_props dictionary.

"""
from __future__ import absolute_import
from __future__ import print_function
import collections
import shelve
import copy
import os
import psi4
import p4util
from p4const import *



def generate_inputs(db,name,displacement_type=1):
    """
        Generates the input files in each sub-directory of the
        distributed finite differences property calculation.

    name: ( string ) method name passed to calling driver,
    db:   (database) The database object associated with this property
          calculation. On exit this db['inputs_generated'] has been set True

    Returns: nothing
    Throws: Exception if the number of atomic displacements is not correct.
    """
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()
    eq_geom = molecule.geometry()
    displacement_geoms = psi4.atomic_displacements(molecule)

    if displacement_type == 2:
        additional_geoms = psi4.mixed_atomic_displacements(molecule)
        for x in additional_geoms:
            displacement_geoms.append(x)

    # Sanity Check
    # there should be 3 cords * natoms *2 directions (+/-)
    # if not (6 * natom) == len(displacement_geoms):
    #     raise p4util.ValidationError('The number of atomic displacements should be 6 times'
    #                     ' the number of atoms!')

    displacement_names = db['job_status'].keys()

    for n, entry in enumerate(displacement_names):
        if not os.path.exists(entry):
            os.makedirs(entry)

        # Setup up input file string
        inp_template = 'molecule {molname}_{disp}'
        inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
        molecule.set_geometry(displacement_geoms[n])
        molecule.fix_orientation(True)
        molecule.fix_com(True)
        inputfile = open('{0}/input.dat'.format(entry), 'w')
        inputfile.write("# This is a psi4 input file auto-generated for"
            "computing properties by finite differences.\n\n")
        inputfile.write(
            inp_template.format(
                molname=molecule.name(),
                disp=entry,
                molecule_info=molecule.create_psi4_string_from_molecule(),
                options=p4util.format_options_for_input(),
                jobspec=db['prop_cmd']))
        inputfile.close()

    if db['do_eq_point']:
        if not os.path.exists('eq'):
            os.makedirs('eq')
        inp_template = 'molecule {molname}_{disp}'
        inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
        molecule.set_geometry(eq_geom)
        molecule.fix_orientation(True)
        molecule.fix_com(True)
        inputfile = open('{0}/input.dat'.format('eq'), 'w')
        inputfile.write("# This is a psi4 input file auto-generated for"
        "computing properties by finite differences.\n\n")
        inputfile.write(
            inp_template.format(
                molname=molecule.name(),
                disp=0,
                molecule_info=molecule.create_psi4_string_from_molecule(),
                options=p4util.format_options_for_input(),
                jobspec=db['prop_cmd']))
        inputfile.close()
    db['inputs_generated'] = True



    # END generate_inputs

def initialize_job_status(database,displacement_type=1):
    """
        Initialize the ordered dict containing job_names and their status.

    database: (database) the databse object passed from the caller
    displacement_type: (int) the order of the partial derivatives being
        computed defaults to 1, which will generate list of singly displaced
        geometries. Currently the only additional accepted value is 2 which
        will generate a list of the singly displaced geoms, and for every geom
        that is singly displaced every coordinate greater than or equal to it
        is also displaced. Coordinates indexed by 0-3*natom.

    Returns: nothing
    Throws: ValidationError if the requested displacement type is unavailable.
    """
    database['job_status'] = collections.OrderedDict()
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()
    coordinates = ['x', 'y', 'z']
    step_direction = ['p', 'm']
    count = 0
    for atom in range(1,natom+1):
        for idx,coord in enumerate(coordinates):
            for step in step_direction:
                job_name_base = "{}_{}_{}".format(atom,coord,step)
                database['job_status'].update(
                    {job_name_base: 'not_started'})
                print( "added element {}: {}".format(count,job_name_base) )
                count +=1
    if displacement_type == 2:
        for i in range(0,natom*3):
            coord1 = i%3
            atom1 = i/3
            for step1 in step_direction:
                for j in range(0,i):
                    atom2 = j/3
                    coord2 = j%3
                    for step2 in step_direction:
                        job_name = "{}_{}_{}_{}_{}_{}".format(
                        atom1+1,coordinates[coord1],step1,
                        atom2+1,coordinates[coord2],step2)
                        print( "added element {}: {}".format(count,job_name) )
                        database['job_status'].update(
                        {job_name: 'not_started'})
                        count +=1


def initialize_database(database, name, prop, properties_array,
        additional_kwargs=None, displacement_type =1):
    """
        Initialize the database for computation of some property
        using distributed finite differences driver

    database: (database) the database object passed from the caller
    name:  (string) name as passed to calling driver
    prop: (string) the property being computed, used to add xxx_computed flag
        to database
    prop_array: (list of strings) properties to go in
        properties kwarg of the property() cmd in each sub-dir
    additional_kwargs: (list of strings) *optional*
        any additional kwargs that should go in the call to the
        property() driver method in each subdir
    displacement_type: (int) type of displacements needed depends on the
        property being computed, and therefore is chosen by the programmer of the
        calling driver, passed to initialize_job_status


    Returns: nothing
    Throws: nothing
    """
    database['inputs_generated'] = False
    database['jobs_complete'] = False
    prop_cmd ="property('{0}',".format(name)
    prop_cmd += "properties=[ '{}' ".format(properties_array[0])
    if len(properties_array) > 1:
        for element in properties_array[1:]:
            prop_cmd += ",'{}'".format(element)
    prop_cmd += "]"
    if additional_kwargs is not None:
        for arg in additional_kwargs:
            prop_cmd += ", {}".format(arg)
    prop_cmd += ")"
    database['prop_cmd'] = prop_cmd
    # Populate the job_status dict
    if displacement_type == 2:
        database['do_eq_point'] = True
        database['eq_point_status'] = 'not_started'
        database['eq_point_{}'.format(prop)] = 0.00
    else:
        database['do_eq_point'] = False
    initialize_job_status(database, displacement_type)

    database['{}_computed'.format(prop)] = False

    # END initialize_database()


def stat(db):
    """
        Checks displacement sub_directories for the status of each
        displacement computation

    db: (database) the database storing information for this distributed
        property calculation

    Returns: nothing
    Throws: nothing
    """
    n_finished = 0
    for job, status in db['job_status'].items():
        if status == 'finished':
            n_finished += 1
        elif status in ('not_started', 'running'):
            try:
                with open("{}/output.dat".format(job)) as outfile:
                    outfile.seek(-150, 2)
                    for line in outfile:
                        if 'Psi4 exiting successfully' in line:
                            db['job_status'][job] = 'finished'
                            n_finished += 1
                            break
                        else:
                            db['job_status'][job] = 'running'
            except:
                pass
    # check all jobs done?
    if db['do_eq_point']:
        if db['eq_point_status'] in ('not_started','running'):
            try:
                with open("eq/output.dat") as eq_outfile:
                    eq_outfile.seek(-150,2)
                    for line in eq_outfile:
                        if 'Pis4 exiting successfully' in line:
                            db['eq_point_status'] = 'finished'
                            break
                        else:
                            db['eq_point_status'] = 'running'
            except:
                pass
    if n_finished == len(db['job_status'].keys()):
        if db['do_eq_point']:
            if db['eq_point_status'] == 'finished':
                db['jobs_complete'] = True
            else:
                return
        else:
            db['jobs_complete'] = True



    # END stat()
