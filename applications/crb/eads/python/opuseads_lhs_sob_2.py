# -*- mode: python; coding: utf-8 -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#        Date: 2010-03-15
#
#   Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
# \file heat1dfemtest.py
# \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
# \author Otmane Souhar <otmane.souhar@imag.fr>
# \date 2010-03-15
#

import os
from openturns import *
from math import *
from pylab import *
from matplotlib.numerix import *
import random
import csv

# load monolithic fem implementation
opuseadspfem = NumericalMathFunction("opuseadspfem")

S1_pfem=[]
S2_pfem=[]

A=[]

f_param = open('plhs.dat', 'r')

s_pfem_Kic_D  = open("output_plhs_pfem_Kic_D_py.txt", "w")
s_pfem_Kic_r  = open("output_plhs_pfem_Kic_r_py.txt", "w")
s_pfem_Kic_ea = open("output_plhs_pfem_Kic_ea_py.txt", "w")
s_pfem_D_r    = open("output_plhs_pfem_D_r_py.txt", "w")
s_pfem_D_ea   = open("output_plhs_pfem_D_ea_py.txt", "w")
s_pfem_r_ea   = open("output_plhs_pfem_r_ea_py.txt", "w")
print "-------------------------------------------\n"

print "CALL OF PFEMTEST"

print "-------------------------------------------\n"
inP = NumericalPoint(opuseadspfem.getInputDimension())
#inP[0] = 2   # kIC : thermal conductivity (default: 2)
#inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
#inP[2] = 1e6  # Q : heat flux (default: 1e6)
#inP[3] = 100  # r : conductance (default: 100)
#inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
#inP[5] = 1 # meshsize times 1e-3 (default: 1)
#inP[6] = 2 # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
#print "-------------------------------------------\n"

N=1000

for line in f_param:
        A.append(line.split())
f_param.close()
A = array(A)

for i in range(0,N):
        inP[0] = 2.0e+0 
        inP[1] = 7.0e-3 
        inP[2] = 1.0e+6 
        inP[3] = float(A[i,2])
        inP[4] = float(A[i,3])
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_Kic_D.write(s1_p)
        s_pfem_Kic_D.write(' ')
        s_pfem_Kic_D.write(s2_p)
        s_pfem_Kic_D.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_Kic_D.close()

print "-------------------------------------------\n"


for i in range(0,N):
        inP[0] = 2.0e+0 
        inP[1] = float(A[i,1])
        inP[2] = 1.0e+6 
        inP[3] = 1.0e+2 
        inP[4] = float(A[i,3])
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_Kic_r.write(s1_p)
        s_pfem_Kic_r.write(' ')
        s_pfem_Kic_r.write(s2_p)
        s_pfem_Kic_r.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_Kic_r.close()

print "-------------------------------------------\n"
for i in range(0,N):
        inP[0] = 2.0e+0 
        inP[1] = float(A[i,1])
        inP[2] = 1.0e+6  
        inP[3] = float(A[i,2])
        inP[4] = 4.0e-3 
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_Kic_ea.write(s1_p)
        s_pfem_Kic_ea.write(' ')
        s_pfem_Kic_ea.write(s2_p)
        s_pfem_Kic_ea.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_Kic_ea.close()

print "-------------------------------------------\n"
for i in range(0,N):
        inP[0] = float(A[i,0])
        inP[1] = 7.0e-3 
        inP[2] = 1.0e+6 
        inP[3] = 1.0e+2 
        inP[4] = float(A[i,3])
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_D_r.write(s1_p)
        s_pfem_D_r.write(' ')
        s_pfem_D_r.write(s2_p)
        s_pfem_D_r.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_D_r.close()

print "-------------------------------------------\n"
for i in range(0,N):
        inP[0] = float(A[i,0])
        inP[1] = 7.0e-3 
        inP[2] = 1.0e+6
        inP[3] = float(A[i,2])
        inP[4] = 4.0e-3
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_D_ea.write(s1_p)
        s_pfem_D_ea.write(' ')
        s_pfem_D_ea.write(s2_p)
        s_pfem_D_ea.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_D_ea.close()

print "-------------------------------------------\n"
for i in range(0,N):
        inP[0] = float(A[i,0])
        inP[1] = float(A[i,1])
        inP[2] = 1.0e+6
        inP[3] = 1.0e+2 
        inP[4] = 4.0e-3 
        inP[5] = 1      # meshsize times 1e-3 (default: 1)
        inP[6] = 2      # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        s1_p=str(outP_pfem[0])
        s2_p=str(outP_pfem[1])
        s_pfem_r_ea.write(s1_p)
        s_pfem_r_ea.write(' ')
        s_pfem_r_ea.write(s2_p)
        s_pfem_r_ea.write('\n')
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))

print "S1_pfem = ", S1_pfem
print "S2_pfem = ", S2_pfem

s_pfem_r_ea.close()

print "-------------------------------------------\n"
