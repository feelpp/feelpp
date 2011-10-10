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

# load monolithic fem implementation
opuseadscrb = NumericalMathFunction("opuseadscrb")


print "-------------------------------------------\n"
inP = NumericalPoint(opuseadscrb.getInputDimension())
inP[0] = 10   # kIC : thermal conductivity (default: 2)
inP[1] = 7e-3 # D : fluid flow rate (default: 5e-3)
inP[2] = 1e6  # Q : heat flux (default: 1e6)
inP[3] = 1e6  # r : conductance (default: 100)
inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)
inP[5] = 1 # meshsize times 1e-3 (default: 1)
inP[6] = 100 # Nwn
inP[7] = 0.01 # maxerror
inP[8] = 3 # errortype

print "-------------------------------------------\n"
N=1
for i in range(0,N):
#inP[0]=exp(log(0.2)+i*log(150)/(N-1))
	inP[0]=10
	print "mu=",inP
	outP = opuseadscrb(inP)
	print "outP = ", outP # s1, s2
print "-------------------------------------------\n"




# Local Variables:
# indent-tabs-mode: t
# End:


