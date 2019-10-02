# -*- Python -*-
#
#  @file  laplacian1dP1.py
#  @brief A test file for the wrapper code
#
#  (C) Copyright 2010 Universite Joseph Fourier
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  This script imports the OpenTURNS environment in Python, loads this wrapper
#  and the wrapped code, and then calls it.
#
import os
from openturns import *
from math import *

# load monolithic fem implementation
laplacian1dP1 = NumericalMathFunction("laplacian1dP1")
# print "laplacian1dP1 = ", laplacian1dP1

print "-------------------------------------------\n"
inP = NumericalPoint(laplacian1dP1.getInputNumericalPointDimension())
inP[0] = 0.1   # a


print "-------------------------------------------\n"
outP = laplacian1dP1(inP)
print "outP = ", outP # s
print "-------------------------------------------\n"

