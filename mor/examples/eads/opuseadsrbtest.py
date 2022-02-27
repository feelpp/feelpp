# -*- Python -*-
#
#  @file  test_opuseadsrb.py
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
from openturns import *

f = NumericalMathFunction("opuseadsrb")
print "f = ", f

inP = NumericalPoint(f.getInputNumericalPointDimension())
inP[0] = 10 # kIC
inP[1] = 5e-3 # D
print "inP = ", inP

outP = f(inP)
print "outP = ", outP # s1, s2
