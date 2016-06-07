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
Module of helper functions for distributed ccresponse computations.

Defines functions for retrieving data computed at displaced geometries.
"""
from __future__ import absolute_import
from __future__ import print_function
import collections
import shelve
import copy
import os
import psi4
import p4util
import numpy as np


def collect_displaced_matrix_data(db, signature, row_dim):
    """
        Gathers a list of tensors, one at each displaced geometry.

    db: (database) the database object for this property calculation
    signature: (string) The string that notifies the matrix reader that the
        targeted tensor data begins.
    row_dim: the expected number of rows that this value should be printed
        across in the file

    Returns a 2d list result[i][j]:
        i: indexes displacements
        j: indexes elements of the flattened tensor at some displacement
    Throws: none
    """
    result = []
    for job in db['job_status']:
        with open('{}/output.dat'.format(job)) as outfile:
            result.append(parse_geometry_matrix_data(outfile, signature, row_dim))

    return result

    # END collect_displaced_matrix_data()


def parse_geometry_matrix_data(outfile, matrix_name, row_tot):
    """
        Parses data from a 3 by n  matrix printed to a file

    outfile: ( file ) handle open in read mode, where the data should be found
    matrix_name: ( string ) that indicates the matrix data is found on the lines
        below
    row_tot: ( int ) indicates the number of lines that the matrix data should
        be printed across in the file

    Returns: matrix_data a list of matrix elements, len = 3*row_tot

    Throws: ParsingError (Collecting matrix data failed) if
            It can't find matrix_header in the file.
            It found matrix_header, but no data.
            It found matrix_header, and data but the number of elements is
            incorrect.

    """
    collect_matrix = False
    n_rows = 0
    n_tries = 0
    matrix_data = []
    for line in outfile:
        if matrix_name in line:
            collect_matrix = True
        if collect_matrix and (n_rows < row_tot):
            try:
                n_tries += 1
                if n_tries > (row_tot + 13):
                    raise ParsingError('{} Matrix was unreadable. Scanned {}'
                                    'lines.'.format(matrix_name, n_tries))
                else:
                    (index, x, y, z) = line.split()
                    matrix_data.append(float(x))
                    matrix_data.append(float(y))
                    matrix_data.append(float(z))
                    n_rows += 1
            except:
                pass
        if (n_rows == row_tot) and (len(matrix_data) != 3 * row_tot):
            raise p4util.ParsingError('Collecting {} data failed!'
                            '\nExpected {} elements but only captured {}'.format(
                                matrix_name, 3 * row_tot, len(matrix_data)))
        if len(matrix_data) == 3 * row_tot:
            return matrix_data

    raise p4util.ParsingError('data for {}  was not found in the output file, '
                    'but it was marked for collection. Check output files '
                    'in displacement sub-dirs!'.format(matrix_name))

    # END parse_geometry_matrix_data()

def parse_hessian_matrix(db, signature, row_dim, natom):
    """
        Searches for file15.dat in the working directory. If found parses the
        data line by line into a 2d array (3natom by 3natom). And returns to
        the caller.
    """
    hessian_read_data = []
    mol = psi4.get_active_molecule()
    try:
        with open('file15.dat') as hessF:
            for line in hessf:
                string_ln = line.strip().split()
                data_ln = [float(x) for x in string_ln]
                hessain_read_data.append(data_ln)

    except IOError:
        # python 2 raises this type
        raise p4util.ParsingError(
                "The hessian data {} could not be found be"
                " sure the hessian has been generated!\n".format('file15.dat'))
    except FileNotFoundError:
        # python 3 raises this type
        raise p4util.ParsingError(
                "The hessian data {} could not be found be"
                " sure the hessian has been generated!\n".format('file15.dat'))

    except:
        raise p4util.ParsingError(
            "A problem occurred while reading the hessian data!\n")

    natom = mol.natom()
    hessian = np.zeros((3*natom,3*natom))
    for i in range(0,3*natom):
        for j in range(0,natom):
            data = hessian_read_data.pop(0)
            hessian[i][3*j+1] = data[0]
            hessian[i][3*j+1] = data[1]
            hessian[i][3*j+1] = data[2]

    return hessian






    # step = psi4.get_local_option('FINDIF', 'DISP_SIZE')

    # alpha_list = []
    # for job in db['job_status']:
    #     with open('{}/output.dat'.format(job)) as outfile:
    #         tensor = parse_geometry_matrix_data(outfile, signature, row_dim)
    #         alpha_list.append(sum( [tensor[3*x] for x in xrange(row_dim)] ))

    # # For the time, lets assume that alpha_list contains alphas in the order:
    # # i) Two displacements for same homogenous derivatives; and then
    # # ii) Four displacements for heterogenous ones after
    # # iii) Equilibrium displacement in the end

    # eq_alpha = alpha_list.pop(-1)

    # result = numpy.zeros((3*natom, 3*natom))
    # for i in xrange(3*natom):
    #     result[i][i] = (alpha_list[2*i] - 2*eq_alpha + alpha_list[2*i + 1])/step/step
    # index = 6*natom
    # for i in xrange(3*natom):
    #     for j in xrange(i):
    #         result[i][j] = (alpha_list[index] - alpha_list[index+1]
    #                                     - alpha_list[index+2] + alpha_list[index+3])/4.0/step/step
    #         result[j][i] = result[i][j]

    # return result
