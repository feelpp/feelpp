# -*- mode: python -*-
#
#  This file is part of the Feel++ library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#        Date: 2012-08-31
#
#   Copyright (C) 2012 Universite de Strasbourg
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
# \file laplacian.py
# \author Christophe Prud'homme <christophe.prudhomme@unistra.fr>
# \date 2012-08-31
#

import os
from openturns import *
from math import *

# load monolithic fem implementation
laplacian = NumericalMathFunction("feelpp_doc_laplacian_ot")


print "-------------------------------------------\n"
inP = NumericalPoint(laplacian.getInputDimension())
inP[0] = 0.2   # hsize: characteristic mesh size


print "-------------------------------------------\n"
print "mu=",inP
outP = laplacian(inP)
print "outP = ", outP # s1, s2
print "-------------------------------------------\n"





# Local Variables:
# indent-tabs-mode: t
# End:


