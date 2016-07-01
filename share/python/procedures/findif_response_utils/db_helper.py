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
from .matrix_from_file_helper import file15_matrix

def generate_inputs(db, database_kwargs=None):
    """
        Generates the input files in each sub-directory of the
        distributed finite differences property calculation.

    name: ( string ) method name passed to calling driver,
    db:   (database) The database object associated with this property
          calculation. On exit this db['inputs_generated'] has been set True

    Returns: nothing
    """
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()
    eq_geom = molecule.geometry()
    dir_name = db['dir_name'];
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    single_displacement_geoms = []
    if database_kwargs['disp_mode']=='cartesian':
        single_displacement_geoms = psi4.atomic_displacements(molecule)
    elif database_kwargs['disp_mode']=='normal':
        hessmat = file15_matrix()
        disp_w_size = psi4.normal_mode_rms_amp_displacements(molecule,hessmat)
        for (disp, size) in disp_w_size:
            vec_m = vec_p = disp.clone()
            vec_p.scale(size)
            vec_m.scale(-1.0*size)
            vec_p.add(eq_geom)
            vec_m.add(eq_geom)
            single_displacement_geoms.extend((vec_p, vec_m))
    else:
        pass

    single_displacement_names = db['jobs']['single_displacements']['job_status'].keys()

    if 'single_displacements' in db['jobs'].keys():
        for n, entry in enumerate(single_displacement_names):
            if not os.path.exists(os.path.join(os.path.join(os.getcwd(), dir_name), entry)):
                os.makedirs(os.path.join(os.path.join(os.getcwd(), dir_name), entry))

            # Setup up input file string
            inp_template = 'molecule {molname}_{disp}'
            inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
            molecule.set_geometry(single_displacement_geoms[n])
            molecule.fix_orientation(True)
            molecule.fix_com(True)
            inputfile = open('{0}/{1}/input.dat'.format(dir_name, entry), 'w')
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

    if 'mixed_displacements' in db['jobs'].keys():

        mixed_displacement_names = db['jobs']['mixed_displacements']['job_status'].keys()
        mixed_displacement_geoms = []
        if database_kwargs['disp_mode']=='cartesian':
            mixed_displacement_geoms = psi4.mixed_atomic_displacements(molecule)
        elif database_kwargs['disp_mode']=='normal':
            raise AssertionError("Mixed displacements in normal mode not yet supported")
        else:
            pass

        for n, entry in enumerate(mixed_displacement_names):
            if not os.path.exists(os.path.join(os.path.join(os.getcwd(), dir_name), entry)):
                os.makedirs(os.path.join(os.path.join(os.getcwd(), dir_name), entry))

            inp_template = 'molecule {molname}_{disp}'
            inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
            molecule.set_geometry(mixed_displacement_geoms[n])
            molecule.fix_orientation(True)
            molecule.fix_com(True)
            inputfile = open('{0}/{1}/input.dat'.format(dir_name, entry), 'w')
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

    if 'eq_point' in db['jobs'].keys():
        if not os.path.exists(os.path.join(os.path.join(os.getcwd(), dir_name), 'eq')):
            os.makedirs(os.path.join(os.path.join(os.getcwd(), dir_name), 'eq'))
        inp_template = 'molecule {molname}_{disp}'
        inp_template += ' {{\n{molecule_info}\n}}\n{options}\n{jobspec}\n'
        molecule.set_geometry(eq_geom)
        molecule.fix_orientation(True)
        molecule.fix_com(True)
        inputfile = open('{0}/{1}/input.dat'.format(dir_name, 'eq'), 'w')
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

def initialize_job_status(database,database_kwargs=None):
    """
    Initialize the ordered dict containing job_names and their status.

    database: (database) the databse object passed from the caller

    Returns: nothing
    """
    molecule = psi4.get_active_molecule()
    natom = molecule.natom()
    # create job_status dict for each job type
    for jb_type in database['jobs']:
        database['jobs'][jb_type].update(
            {'job_status': collections.OrderedDict()}
            )

    # The ordering of the job names here also matches up with the
    # return ordering from the functions which create the displaced
    # geometries. I would prefer those lookups created a dict (std::map)
    # but for now that would be too much to overhaul
    # always have single displacements so set them up

    if database_kwargs['disp_mode'] == 'cartesian':
        coordinates = ['x', 'y', 'z']
        step_direction = ['p', 'm']

        for atom in range(1,natom+1):
            for idx,coord in enumerate(coordinates):
                for step in step_direction:
                    job_name= "{}_{}_{}".format(atom,coord,step)
                    database['jobs']['single_displacements']['job_status'].update(
                        {job_name: 'not_started'})

        # displacements generated in the same order (convenient later when writing
        # input files!). Each Cartesian coordinate coord1 (3* natom total) are displaced
        # then for every coordinate from 0 to coord1 inclusive are displaced. This
        # generated all displacements needed to compute the lower triangle of the
        # 2nd derivatives, less the main diagonal. The main diagonal is computed
        # using the singly displaced coordinates in the single_displacements list,
        # and the eq point.
        if 'mixed_displacements' in database['jobs'].keys():
            for i in range(0,natom*3):
                coord1 = i%3
                atom1 = i/3
                for step1 in step_direction:
                    # we don't need cases where i and j are equal so this is fine
                    # the case where they are equal would be for the diagonal
                    # elements of the 4th derivatives using the form we are here
                    # these points may be needed if we change to a different scheme
                    for j in range(0,i):
                        atom2 = j/3
                        coord2 = j%3
                        for step2 in step_direction:
                            job_name = "{}_{}_{}_{}_{}_{}".format(
                            atom1+1,coordinates[coord1],step1,
                            atom2+1,coordinates[coord2],step2)
                            database['jobs']['mixed_displacements']['job_status'].update(
                            {job_name: 'not_started'})


        if 'eq_point' in database['jobs'].keys():
            database['jobs']['eq_point']['job_status'].update(
                {'eq': 'not_started'}
                )

    elif database_kwargs['disp_mode'] == 'normal':
        num_modes = 3*natom - 6
        # modes = range(1,num_modes + 1)
        step_direction = ['p', 'm']

        for curr_mode in range(1, num_modes + 1):
            for step in step_direction:
                job_name= "mode{}_{}".format(curr_mode, step)
                database['jobs']['single_displacements']['job_status'].update(
                    {job_name: 'not_started'}
                )

        # Same way as before. One could require mixed displacements in normal modes,
        # sometime in future.
        if 'mixed_displacements' in database['jobs'].keys():
            raise AssertionError('Mixed displacement not yet supported for normal modes')
        else:
            pass

        if 'eq_point' in database['jobs'].keys():
            database['jobs']['eq_point']['job_status'].update(
                {'eq': 'not_started'}
            )

    else:
        pass


def initialize_database(database, name, prop, properties_array, dir_name,
        database_kwargs = None):
    """
        Initialize the database for computation of some property
        using distributed finite differences driver

    database: (database) the database object passed from the caller
    name:  (string) name as passed to calling driver
    prop: (string) the property being computed, used to add xxx_computed flag
        to database
    prop_array: (list of strings) properties to go in
        properties kwarg of the property() cmd in each sub-dir
    database_kwargs: (list of strings) *optional*
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
    database['dir_name'] = dir_name
    prop_cmd ="property('{0}',".format(name)
    prop_cmd += "properties=[ '{}' ".format(properties_array[0])
    if len(properties_array) > 1:
        for element in properties_array[1:]:
            prop_cmd += ",'{}'".format(element)
    prop_cmd += "]"
    # if database_kwargs is not None:
    #     for arg in database_kwargs:
    #         prop_cmd += ", {}".format(arg)
    prop_cmd += ")"
    database['prop_cmd'] = prop_cmd

    # Populate the jobs dict, create dict for each 'type' of job
    # so far this organization makes sense for all routines
    # this may need to change later
    database['jobs'] = collections.OrderedDict()

    if 'single' in database_kwargs['disp_type']: 
        database['jobs'].update(
            {'single_displacements': collections.OrderedDict()}
            )
    # if we need 2nd derivatives, generate their dict
    # also one for the eq point
    if 'mixed' in database_kwargs['disp_type']: 
        database['jobs'].update(
            {'mixed_displacements': collections.OrderedDict()}
            )

    if 'eq' in database_kwargs['disp_type']:
        database['jobs'].update(
            {'eq_point': collections.OrderedDict()}
            )

    initialize_job_status(database, database_kwargs)

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
    n_total = 0
    dir_name = db['dir_name']
    for job_type in db['jobs'].keys():
        n_total += len(db['jobs'][job_type]['job_status'])
        for job, status in db['jobs'][job_type]['job_status'].items():
            if status == 'finished':
                n_finished += 1
            elif status in ('not_started', 'running'):
                try:
                    with open("{0}/{1}/output.dat".format(dir_name, job)) as outfile:
                        outfile.seek(-150, 2)
                        for line in outfile:
                            if 'Psi4 exiting successfully.' in line:
                                db['jobs'][job_type]['job_status'][job] = 'finished'
                                n_finished += 1
                                break
                            else:
                                db['jobs'][job_type]['job_status'][job] = 'running'
                except:
                    pass
    if n_finished == n_total:
        db['jobs_complete'] = True

    return (n_finished,n_total)


    # END stat()

# ################################
# ###                          ###
# ###    DATABASE STRUCTURE    ###
# ###                          ###
# ################################
# Dict of dicts
# 'inputs_generated' ==> (boolean)
# 'jobs_complete'    ==> (boolean)
# 'prop_cmd'         ==> (string, "property(name, [properties], addn_props)"
# 'dir_name'         ==> <string>
# 'jobs'             ==> (OrderedDict)
#       'single_displacements' :: OrderedDict
#           'job_status' :: OrderedDict (keys 1_x_p,... value 'not_started',...)
#       'mixed_displacements' :: OrderedDict (if type=2)
#           'job_status' :: OrderedDict (keys 1_x_p,... value 'not_started',...)
#       'eq_point'            :: OrderedDict (if type=2)
#           'job_status' :: OrderedDict (keys 1_x_p,... value 'not_started',...)
# '<prop>_computed'  ==> (boolean)
#
##################################
