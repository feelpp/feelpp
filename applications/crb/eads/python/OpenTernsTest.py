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

################################
###   Function 'opus_mfem'   ###
################################
# Create here the python lines to define the implementation of the function

# In order to be able to use that function with the openturns library, it is necessary to define aclass which derives from OpenTURNSPythonFunction

class modelePYTHON(OpenTURNSPythonFunction):
# That following method defines the input size (5) and the output size(7): This implementation is for the parametrized FEM (PFEM)
        def __init__(self):
          OpenTURNSPythonFunction.__init__(self, 5, 7)
# That following method gives the implementation of modelePYTHON
        def f(self,X):
	  v = []
	  for i in range(len(X)):
	    v.append(X[i])
	  v.append(1.) # meshsize times 1e-3 (default: 1)
	  v.append(2) # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	  return v

# Use that function defined in the script python with the openturns library
# Create a NumericalMathFunction from modelePYTHON
 
# Load monolithic fem implementation
opuseadspfem = NumericalMathFunction("opuseadspfem")
eads_mfem = NumericalMathFunction("opuseadsmfem")


