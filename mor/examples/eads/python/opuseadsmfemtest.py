# -*- Python -*-
#
#  @file  test_opuseadsfem.py
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
eads_mfem = NumericalMathFunction("opuseadsmfem")
# print "f = ", f

# base_dir=os.environ["HOME"]+"/life/opus/D_"+str(inP[1])+"/kIC_"+str(inP[0])+"/"
# base_dir=os.environ["HOME"]+"/life/opus/"
#os.makedirs(base_dir)
# os.chdir(base_dir)

print "-------------------------------------------\n"
inP = NumericalPoint(eads_mfem.getInputDimension())
inP[0] = 10   # kIC : thermal conductivity (default: 2)
inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
inP[2] = 1e6  # Q : heat flux (default: 1e6)
inP[3] = 100  # r : conductance (default: 100)
inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
inP[5] = 1 # meshsize times 1e-3 (default: 1)
inP[6] = 2 # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)


print "-------------------------------------------\n"
N=10
for i in range(0,N):
    inP[0]=exp(log(0.2)+i*log(150)/(N-1))
    print "mu=",inP
    outP = eads_mfem(inP)
    print "outP = ", outP # s1, s2
print "-------------------------------------------\n"


