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
from __future__ import absolute_import
from __future__ import print_function

from p4const import *
import p4util
import psi4
import collections
import shelve
import copy
import sys
import inspect
import os
# Relative hack for now
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0], "../")))
sys.path.append(path_dir)
from . import findif_response_utils


# maybe rename this rotation_findif_correction ??
def rotation_findif_correction(name, **kwargs):
    """
        Main driver for managing Zero Point Correction to Optical activity
        computations with CC response theory.

        Uses distributed finite differences approach -->
            1. Sets up a database to keep track of running/finished/waiting
                computations.
            2. Generates separate input files for displaced geometries.
            3. When all displacements are run, collects the necessary information
                from each displaced computation, and computes final result.
    """
# Relevant kwargs will be
    # correction_type: vibave, zpvc
    # displacement_coords: cartesian, normal, normal_rms_amp

    # for now lets assume correction_type is vibave, deal with zpvc later

    # Get list of omega values -> Make sure we only have one wavelength
    # Catch this now before any real work gets done
    omega = psi4.get_option('CCRESPONSE', 'OMEGA')
    if len(omega) > 2:
        raise p4util.ValidationError('rotation_findif_correction can only be performed for one wavelength.')
    else:
        pass

    ## assume in NM
    nu_hartrees = psi_h*psi_c/(omega[0]*(1e-9))/psi_hartree2J

    mol = psi4.get_active_molecule()

    ## by default use cartesian coords
    if not('disp_mode' in kwargs): 
        kwargs.update({'disp_mode' : 'cartesian'})

    ############################
    ## if one passes disp_mode = 'normal', then the regex(^(no|false|off|0))
    ## in kwargs_lower(in procutil.py) matches the string 'normal' and is 
    ## in turn converted to False.
    ## Work around for now, let us fix that later
    if kwargs['disp_mode'] == False:
        kwargs.update({'disp_mode' : 'normal'})
    ############################

    psi4.print_out(
        'Running rotation_findif_correction computation. Subdirectories for each '
        'required displaced geometry have been created.\n\n')

    # Initialize database
    db = shelve.open('database', writeback=True)
    dbno = 0

    # Check if final result is in here
    # ->if we have already computed property, back up the dict
    # ->copy it setting this flag to false and continue
    # dirs are created in order dir0, dir1...
    # respective daatabases are _bak0, _bak1....
    # current working database is plane database which corresponds to most fresh dir
    # if ('zpvc_computed' in db) and ( db['zpvc_computed'] ):
    #     # increment dbno until backup database exists
    #     while True:
    #         db_bak = shelve.open('.database.bak{}'.format(dbno), writeback=True)
    #         if not ('zpvc_computed' in db_bak):
    #             break
    #         else:
    #             db_bak.close()
    #             dbno += 1

    #     for key,value in db.iteritems():
    #         db_bak[key]=value

    #     db_bak.close()
    #     db['zpvc_computed'] = False
    
    #     dbno += 1
    # else:
    #     db['zpvc_computed'] = False
    # Assertion: dbno = number of dirs created before this

    database_kwargs = {}
    if kwargs['disp_mode'] == 'normal':
        database_kwargs = {'disp_mode' : 'normal', 
                          'disp_type' : ['single', 'eq']
                        }
    elif kwargs['disp_mode'] == 'cartesian':
        database_kwargs = {'disp_mode' : 'cartesian', 
                          'disp_type' : ['single', 'mixed','eq']
                        }

    # dbno = number of dirs already made
    if 'inputs_generated' not in db:
        findif_response_utils.initialize_database(db,name,"zpvc_rotation",
                ["rotation"],"dir{}".format(dbno), database_kwargs)

    # Generate input files
    if not db['inputs_generated']:
        findif_response_utils.generate_inputs(db, database_kwargs)

    # Check job status
    if db['inputs_generated'] and not db['jobs_complete']:
        print('Checking status')
        (n_finished, n_total) = findif_response_utils.stat(db)
        for job_type in db['jobs']:
            print ("Checking {} jobs".format(job_type))
            for job,status in db['jobs'][job_type]['job_status'].items():
                print("{} --> {}".format(job, status))
        print ("{:.2f}% jobs({:d}/{:d}) completed.".format(
            (n_finished*100.0)/n_total, n_finished, n_total)
        )

    # Compute ZPVC_rotation
    if db['jobs_complete']:

        mygauge = psi4.get_option('CCRESPONSE', 'GAUGE')
        consider_gauge = {
            'LENGTH': ['Length Gauge'],
            'VELOCITY': ['Modified Velocity Gauge'],
            'BOTH': ['Length Gauge', 'Modified Velocity Gauge']
        }
        gauge_list = ["{} Results".format(x) for x in consider_gauge[mygauge]]

        opt_rot_single = []
        opt_rot_eq = []
        if kwargs['disp_mode'] == 'cartesian':
            opt_rot_mixed = []

        for gauge in consider_gauge[mygauge]:
            # Gather data from single_displacements
            opt_rot_single.append(
                findif_response_utils.collect_displaced_matrix_data(
                    db, 'single_displacements',"Optical Rotation Tensor ({})".format(gauge),
                    3)
                )
            # Gather data from mixed displacements, needed only in case cartesian displacements
            if kwargs['disp_mode'] == 'cartesian':
                opt_rot_mixed.append(
                    findif_response_utils.collect_displaced_matrix_data(
                        db, 'mixed_displacements',"Optical Rotation Tensor ({})".format(gauge),
                        3)
                    )
            # Gather data from equilibrium point
            opt_rot_eq.append(
                findif_response_utils.collect_displaced_matrix_data(
                    db, 'eq_point',"Optical Rotation Tensor ({})".format(gauge),
                    3)
                )

        # Structure of opt_rot_single:
        #     {list of guages}->{list of tensors}->{list of 9 floats}
        alpha_single = [
            [
                sum([float(tensor[4*index]) for index in range(3)]) for tensor in guage
            ]
        for guage in opt_rot_single
        ]

        if kwargs['disp_mode'] == 'cartesian':
            alpha_mixed = [
                [
                    sum([float(tensor[4*index]) for index in range(3)]) for tensor in guage
                ]
            for guage in opt_rot_mixed
            ]

        alpha_eq = [
            [
                sum([float(tensor[4*index]) for index in range(3)]) for tensor in guage
            ]
        for guage in opt_rot_eq
        ]


        for i,curr_guage in enumerate(consider_gauge[mygauge]):

            if kwargs['disp_mode'] == 'cartesian':
                D_2 = []
                step = psi4.get_local_option('FINDIF', 'DISP_SIZE')
                
                for j in range(len(opt_rot_single)/2):
                # j enumerates the 3n coordinates
                # '2j' is displacement in +ve direction in coordinate 'j'
                # '2j+1' in the -ve direction
                    curr_row = []
                    for k in range(j+1):
                        val = 0.0

                        if (k == j):
                            numr = alpha_single[i][2*j]
                            numr += alpha_single[i][2*j + 1]
                            numr -= 2*alpha_eq[i][0]

                            val = numr/(step*step)
                        else:
                            numr = alpha_mixed[i][idx(2*j, 2*k)]
                            numr += alpha_mixed[i][idx(2*j + 1, 2*k + 1)]
                            numr -= alpha_mixed[i][idx(2*j, 2*k + 1)]
                            numr -= alpha_mixed[i][idx(2*j + 1, 2*k)]

                            val = numr/(4*step*step)

                        curr_row.append(val)

                    ## D_2 ==> lower triangle of D matrix
                    D_2.append(curr_row)

                ## print(D_2)
                ## Call function and print result
                ## psi4.rotation_vibave_cartesian(mol, D_2)
                ## ideal name ^
                psi4.rotation_vibave_cartesian(mol, D_2)
                ## hack for now ^
                ## Do the necessary parts over here

            if kwargs['disp_mode'] == 'normal':
                deriv_list = []
                step = psi4.get_local_option('FINDIF', 'DISP_SIZE')

                const_pi = 3.14159265358979323846264338327950288

                ## fix this v
                hessmat = findif_response_utils.file15_matrix()
                disp_sizes = psi4.normal_mode_rms_amp_displacements(mol,hessmat)
                disp_amps = [ size*psi_bohr2angstroms for (disp, size) in disp_sizes ]
                ## fix this ^
                ## actually not needed, `delta(x)^2` just cancel each other out
                ## in final expression, so can remove this block

                #calculate factor of conversion to degrees
                bohr2a4 = psi_bohr2angstroms * psi_bohr2angstroms * psi_bohr2angstroms * psi_bohr2angstroms
                m2a = psi_bohr2angstroms * 1.0e-10
                hbar = psi_h/(2.0 * const_pi)
                prefactor = 1.0e-2 * hbar/(psi_c * 2.0 * const_pi * psi_me * m2a * m2a)
                prefactor *= prefactor
                prefactor *= 288.0e-30 * const_pi * const_pi * psi_na * bohr2a4
                print("python prefactor = {}".format(prefactor))
                total_mass = sum([mol.mass(atom_idx) for atom_idx in range(mol.natom())])
                
                print("Total mass {}".format(total_mass))
                if curr_guage == 'Length Gauge':
                    factor = prefactor * (-1.0/3.0/nu_hartrees) * nu_hartrees * nu_hartrees / total_mass
                elif curr_guage == 'Modified Velocity Gauge':
                    factor = prefactor * (-1.0/3.0/nu_hartrees/nu_hartrees) * nu_hartrees * nu_hartrees / total_mass
                print("python factor = {}".format(factor))

                for j in range(len(opt_rot_single[i])/2):
                # j enumerates the 3n coordinates
                # '2j' is displacement in +ve direction in coordinate 'j'
                # '2j+1' in the -ve direction
                    numr = alpha_single[i][2*j]
                    numr += alpha_single[i][2*j+1]
                    numr -= 2*alpha_eq[i][0]

                    val = numr/(disp_amps[j]*step*disp_amps[j]*step)

                    deriv_list.append(val)

                corrections_list = [(0.5)*alpha*delta_x*delta_x for (alpha, delta_x) in zip(deriv_list, disp_amps)]
                correction = sum(corrections_list)
                
                freq_cm = []
                try:
                    freq_cm = [ line.rstrip('\n') for line in open('freq.dat')][3:3+(3*mol.natom())-6]
                    freq_cm = [ float(line.split(':')[1]) for line in freq_cm]
                except:
                    print("file 'freq.dat' not found")

                psi4.print_out("\n\n")
                psi4.print_out("##########################################################\n")
                psi4.print_out("#### {0} Results at {1:6.2f} nm##\n".format(curr_guage, omega[0]))
                psi4.print_out("##########################################################\n")
                psi4.print_out("\n")
                psi4.print_out("---------------------------------------------------------------------------------------------\n")
                psi4.print_out("Mode  freq    alpha(deg/[dm (g/cm^3)])         delta_x           deriv         correction\n")
                psi4.print_out("     (cm-1)   positive        negative        angstrom                      deg/[dm (g/cm^3)]\n")
                psi4.print_out("---------------------------------------------------------------------------------------------\n")
                for idx, size in enumerate(disp_amps):
                    psi4.print_out("{0}\t{1:7.2f}\t{2:12.8f}\t{3:12.8f}\t{4:12.8f}\t{5:13.8f}\t{6:12.8f}\n".format(
                        idx,freq_cm[idx],alpha_single[i][2*idx]*factor,alpha_single[i][2*idx+1]*factor,disp_amps[idx],deriv_list[idx]*factor, corrections_list[idx]*factor)
                    )
                psi4.print_out("---------------------------------------------------------------------------------------------\n")

                psi4.print_out("\n")
                psi4.print_out("Rotation(in deg/[dm (g/cm^3)]):\t{0:12.8f}\n".format(alpha_eq[i][0]*factor))
                psi4.print_out("Correction(in deg/[dm (g/cm^3)]):\t{0:12.8f}\n".format(correction*factor))
                psi4.print_out("##########################################################\n")
                psi4.print_out("\n\n")

        db['zpvc_computed'] = True

    db.close()

# Function that gives index corresponding to (j, k)
# Here j, k correspond to displacements for example, "1_y_m", "1_x_p" etc.
# idx(a, b) = 2((a/2)-1)(a/2) + (a%2)(a-1) + b
def idx(a, b):
    if a < b:
        swap(a,b)
    return (2*((a/2)-1)*(a/2) + (a%2)*(a-1) + b)

#   SAVE this for when multiple wavelengths works
# # Get list of omega values
# omega = psi4.get_option('CCRESPONSE','OMEGA')
# if len(omega) > 1:
#     units = copy.copy(omega[-1])
#     omega.pop()
# else:
#     units = 'atomic'
# wavelength = copy.copy(omega[0])
# # Set up units for scatter.cc
# if units == 'NM':
#     wavelength = (psi_c * psi_h * 1*(10**-9))/(wavelength * psi_hartree2J)
# if units == 'HZ':
#     wavelength = wavelength * psi_h / psi_hartree2J
# if units == 'EV':
#     wavelength = wavelength / psi_hartree2ev
# if units == 'atomic':
#     pass

# view ./findif_response_utils/db_helper.py for database structure
