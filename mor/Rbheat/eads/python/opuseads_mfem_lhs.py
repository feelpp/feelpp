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
from pylab import *
from matplotlib.numerix import *
import random
import csv

# load monolithic fem implementation
eads_mfem = NumericalMathFunction("opuseadsmfem")

S1_mfem=[]
S2_mfem=[]

A=[]

f_param = open('plhs1.dat', 'r')
s_mfem_all = open("output_plhs_mfem_all_py.txt", "w")

print "-------------------------------------------\n"

print "CALL OF MFEMTEST"

print "-------------------------------------------\n"
inP = NumericalPoint(eads_mfem.getInputNumericalPointDimension())
#inP[0] = 2   # kIC : thermal conductivity (default: 2)
#inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
#inP[2] = 1e6  # Q : heat flux (default: 1e6)
#inP[3] = 100  # r : conductance (default: 100)
#inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
#inP[5] = 1 # meshsize times 1e-3 (default: 1)
#inP[6] = 2 # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
#print "-------------------------------------------\n"

N=100

for line in f_param:
        A.append(line.split())
f_param.close()
A = array(A)

print "-------------------------------------------\n"
for i in range(0,N):
        inP[0] = float(A[i,0])
        inP[1] = float(A[i,1])
        inP[2] = 1.0e+6
        inP[3] = float(A[i,2])
        inP[4] = float(A[i,3])
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_mfem = ", inP
	outP_mfem =eads_mfem(inP)
	print "outP_mfem = ", outP_mfem # s1, s2
        s1_m=str(outP_mfem[0])
        s2_m=str(outP_mfem[1])
        s_mfem_all.write(s1_m)
        s_mfem_all.write(' ')
        s_mfem_all.write(s2_m)
        s_mfem_all.write('\n')
        S1_mfem=concatenate((S1_mfem, outP_mfem[0]), axis=None)
        vstack((S1_mfem, S1_mfem))
        S2_mfem=concatenate((S2_mfem, outP_mfem[1]), axis=None)
        vstack((S2_mfem, S2_mfem))

print "S1_mfem = ", S1_mfem
print "S2_mfem = ", S2_mfem

s_mfem_all.close()

print "-------------------------------------------\n"
