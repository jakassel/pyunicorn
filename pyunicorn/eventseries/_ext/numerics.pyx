# -*- coding: utf-8 -*-
#
# This file is part of pyunicorn.
# Copyright (C) 2008--2020 Jonathan F. Donges and pyunicorn authors
# URL: <http://www.pik-potsdam.de/members/donges/software>
# License: BSD (3-clause)
#
# Please acknowledge and cite the use of this software and its authors
# when results are used in publications or published elsewhere.
#
# You can use the following reference:
# J.F. Donges, J. Heitzig, B. Beronov, M. Wiedermann, J. Runge, Q.-Y. Feng,
# L. Tupikina, V. Stolbova, R.V. Donner, N. Marwan, H.A. Dijkstra,
# and J. Kurths, "Unified functional network and nonlinear time series analysis
# for complex systems science: The pyunicorn package"

cimport cython
from cpython cimport bool


import numpy as np
cimport numpy as np

BOOLTYPE = np.uint8
INTTYPE = np.int
INT8TYPE = np.int8
INT16TYPE = np.int16
INT32TYPE = np.int32
FLOATTYPE = np.float
FLOAT32TYPE = np.float32
FLOAT64TYPE = np.float64
ctypedef np.uint8_t BOOLTYPE_t
ctypedef np.int_t INTTYPE_t
ctypedef np.int8_t INT8TYPE_t
ctypedef np.int16_t INT16TYPE_t
ctypedef np.int32_t INT32TYPE_t
ctypedef np.float_t FLOATTYPE_t
ctypedef np.float32_t FLOAT32TYPE_t
ctypedef np.float64_t FLOAT64TYPE_t


#cdef extern from "src_numerics.c":

from libcpp.vector cimport vector

cdef void _calculate_eca():
    param1 = 1
    param2 = 2
    return param1, param2

cdef void _calculate_es():
    param1 = 1
    return param1