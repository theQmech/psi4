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


def run_zpvc_rotation(name, **kwargs):
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

    # Get list of omega values -> Make sure we only have one wavelength
    # Catch this now before any real work gets done
    omega = psi4.get_option('CCRESPONSE', 'OMEGA')
    if len(omega) > 2:
        raise p4util.ValidationError('ZPVC_rotation can only be performed for one wavelength.')
    else:
        pass

    psi4.print_out(
        'Running ZPVC_rotation computation. Subdirectories for each '
        'required displaced geometry have been created.\n\n')

    dbno = 0
    # Initialize database
    db = shelve.open('database', writeback=True)
    # Check if final result is in here
    # ->if we have already computed property, back up the dict
    # ->copy it setting this flag to false and continue
    # if ('zpvc_computed' in db) and ( db['zpvc_computed'] ):
    #     db2 = shelve.open('.database.bak{}'.format(dbno), writeback=True)
    #     dbno += 1
    #     for key,value in db.iteritems():
    #         db2[key]=value

    #     db2.close()
    #     db['zpvc_computed'] = False
    # else:
    #     db['zpvc_computed'] = False

    if 'inputs_generated' not in db:
        findif_response_utils.initialize_database(db,name,"zpvc_rotation",
                ["rotation"],additional_kwargs=None,displacement_type =2 )

    # Generate input files
    if not db['inputs_generated']:
        findif_response_utils.generate_inputs(db,name)

    # Check job status
    if db['inputs_generated'] and not db['jobs_complete']:
        print('Checking status')
        findif_response_utils.stat(db)
        for job_type in db['jobs']:
            print ("Checking {} jobs".format(job_type))
            for job,status in db['jobs'][job_type]['job_status'].items():
                print("{} --> {}".format(job, status))

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
        opt_rot_mixed = []
        opt_rot_eq = []
        for gauge in consider_gauge[mygauge]:
            # Gather data from single_displacements
            opt_rot_single.append(
                findif_response_utils.collect_displaced_matrix_data(
                    db['jobs']['single_displacements']['job_status'].keys(),
                    "Optical Rotation Tensor ({})".format(gauge),
                    3)
                )
            # Gather data from mixed displacements
            opt_rot_mixed.append(
                findif_response_utils.collect_displaced_matrix_data(
                    db['jobs']['mixed_displacements']['job_status'].keys(),
                    ""
                    "Optical Rotation Tensor ({})".format(gauge),
                    3)
                )
            # Gather data from equilibrium point
            opt_rot_eq.append(
                findif_response_utils.collect_displaced_matrix_data(
                    db['jobs']['eq_point']['job_status'].keys(),
                    "Optical Rotation Tensor ({})".format(gauge),
                    3)
                )

        # tmp = [opt_rot_single, opt_rot_mixed, opt_rot_eq]
        # [alpha_single, alpha_mixed, alpha_eq] = 
        #     [ 
        #         [ 
        #             ([(sum([tensor[4*index] for index in xrange(3)]))/3.0 for tensor in guage]) 
        #         for guage in disp_type
        #         ] 
        #     for disp_type in tmp 
        #     ]

# Structure of opt_rot_single:
#     {list of guages}->{list of tensors}->{list of 9 floats}
        alpha_single = [
            [ 
                sum([float(tensor[4*index]) for index in xrange(3)]) for tensor in guage
            ] 
        for guage in opt_rot_single
        ]

        alpha_mixed = [
            [ 
                sum([float(tensor[4*index]) for index in xrange(3)]) for tensor in guage
            ] 
        for guage in opt_rot_mixed
        ]

        alpha_eq = [
            [ 
                sum([float(tensor[4*index]) for index in xrange(3)]) for tensor in guage
            ] 
        for guage in opt_rot_eq
        ]

        step = psi4.get_local_option('FINDIF', 'DISP_SIZE')

        for i in xrange(len(consider_gauge[mygauge])):
            deriv_list = []
            for j in xrange(len(opt_rot_single[i])/2):
            # j enumerates the 3n coordinates
            # '2j' is displacement in +ve direction in coordinate 'j'
            # '2j+1' in the -ve direction
                curr_list = []
                for k in xrange(j+1):
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

                    curr_list.append(val)

                deriv_list.append(curr_list)

            # print(deriv_list)
            ## Call function and print result
            psi4.zpvc_rotation(psi4.get_active_molecule(), deriv_list)
            ## Do the necessary parts over here

    #TODO:
    # - compute 2nd derivatives Cartesian
    # - read hessian
    # - compute normal mode vectors
    # - build rotation-translation projector
    # - transform 2nd derivatives Cartesian to normal modes
    # - project out rotation and translations
    # - compute zpvc to optical rotation
    # - report rotation @ eq-point
    # - report correction alone
    # - report total rotation+correction

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
