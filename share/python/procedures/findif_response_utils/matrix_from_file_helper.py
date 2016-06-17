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

Defines functions that return psi4 matrix types, specific to reading some
particular files
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


def file15_matrix(hessfile='file15.dat'):
    mol = psi4.get_active_molecule()
    natom = mol.natom();
    F = np.zeros((natom*3,natom*3))
    try:
        with open(hessfile) as hessian_file:
            hessian_data = hessian_file.readlines()
            line = 0
            for i in xrange(3*natom):
                for j in xrange(natom):
                    data_line = hessian_data[line]
                    line +=1
                    [F[i,3*j],F[i,(3*j)+1],F[i,(3*j)+2]] = [float(x) for x in
                        data_line.split()]
        return p4util.array_to_matrix(F)
    except IOError:
        # python 2 error
        raise p4util.ParsingError(
                "The hessian data could not be found in {}"
                "be sure the hessian has been generated!\n".format('file15.dat'))

    except FileNotFoundError:
        # python 3 error
        raise p4util.ParsingError(
                "The hessian data could not be found in {}"
                "be sure the hessian has been generated!\n".format('file15.dat'))
