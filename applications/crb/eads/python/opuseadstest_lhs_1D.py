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
# \date 2010-03-15
#

import os
from openturns import *
from math import *
from pylab import *
from matplotlib.numerix import *
import random

# load monolithic fem implementation
opuseadspfem = NumericalMathFunction("opuseadspfem")
eads_mfem = NumericalMathFunction("opuseadsmfem")

D=[]
S1_pfem=[]
S2_pfem=[]
S1_mfem=[]
S2_mfem=[]
Erreur_S1=[]
Erreur_S2=[]

print "-------------------------------------------\n"

print "CALL OF PFEMTEST"

print "-------------------------------------------\n"
inP = NumericalPoint(opuseadspfem.getInputDimension())
inP[0] = 10   # kIC : thermal conductivity (default: 2)
inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
inP[2] = 1e6  # Q : heat flux (default: 1e6)
inP[3] = 100  # r : conductance (default: 100)
inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
#print "-------------------------------------------\n"
N=10
size=1.0/N
minn=0.02
maxx=150
for i in range(0,N):
        minn = i*size
        point = minn + (random.random() * size)
	inP[0]= (point * (maxx - minn)) + minn
        D=concatenate((D, inP[0]), axis=None)
        vstack((D, D)) 
	print "inP_pfem = ", inP
	outP_pfem = opuseadspfem(inP)
	print "outP_pfem = ", outP_pfem # s1, s2
        S1_pfem=concatenate((S1_pfem, outP_pfem[0]), axis=None)
        vstack((S1_pfem, S1_pfem))
        S2_pfem=concatenate((S2_pfem, outP_pfem[1]), axis=None)
        vstack((S2_pfem, S2_pfem))
print "-------------------------------------------\n"

print "CALL OF MFEMTEST"

print "-------------------------------------------\n"
inP = NumericalPoint(eads_mfem.getInputNumericalPointDimension())

inP[0] = 10   # kIC : thermal conductivity (default: 2)
inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
inP[2] = 1e6  # Q : heat flux (default: 1e6)
inP[3] = 100  # r : conductance (default: 100)
inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
inP[5] = 1 # meshsize times 1e-3 (default: 1)
inP[6] = 2 # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
print "-------------------------------------------\n"
for i in range(0,N):
        inP[0]=D[i] 
	#inP[0]=exp(log(0.2)+i*log(150)/(N-1))
	print "inP_mfem = ", inP
        outP_mfem = eads_mfem(inP)
	print "outP_mfem = ", outP_mfem # s1, s2
        S1_mfem=concatenate((S1_mfem, outP_mfem[0]), axis=None)
        vstack((S1_mfem, S1_mfem))
        S2_mfem=concatenate((S2_mfem, outP_mfem[1]), axis=None)
        vstack((S2_mfem, S2_mfem))
print "-------------------------------------------\n"

for i in range(1,N):
         Erreur_S1 = S1_pfem-S1_mfem
         Erreur_S2 = S2_pfem-S2_mfem

print "Erreur_S1 = ", Erreur_S1
print "Erreur_S2 = ", Erreur_S2


plot(D, S1_pfem, 'go-', label='line 1', linewidth=2)

#plot([1,2,3], [1,4,9], 'rs',  label='line 2')
#axis([0, 4, 0, 10])
#legend()

# Local Variables:
# indent-tabs-mode: t
# End:


