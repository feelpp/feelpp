#! /usr/bin/python

# -*- mode: python -*-
#
#  This file is part of the Feel library
#
#  Author(s): Goncalo Pena <gpena@mat.uc.pt>
#        Date: 2010-01-09
#
#   Copyright (C) 2010 University of Coimbra
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
#   University of Coimbra
#
# \file polysetvis.py
# \author Goncalo Pena <gpena@mat.uc.pt>
# \date 2010-01-09
#

# This file draws the equidistributed point set for the several elements and orders


from pyx import *
from scipy import *

import os
import sys
import string

import math
import hypercube
import simplex


from optparse import OptionParser
# Define options
parser = OptionParser()
parser.add_option("--dim", dest="dim", type="int", help="topological dimension of the element", default=2)
parser.add_option("--element_type", dest="element", type="string", help="type of the reference element: simplex or hypercube", default="simplex")
parser.add_option("--order", dest="order", type="int", help="pointset order", default=1)
parser.add_option("--numbering", dest="numbers", type="int", help="adds numbering to points: 0 for no numbering, 1 to add the numbers", default=0)
parser.add_option("--circle_radius", dest="circle_radius", type="float", help="radius of the circles that mark the nodes", default=0.07)
parser.add_option("--size", dest="line_size", type="float", help="length of the side of the element", default=5)
(options, args) = parser.parse_args()


# Verification of the options
if (options.order < 0):
	print "[validation] error : order must be nonnegative"
	sys.exit(-1)

if not (options.numbers == 0 or options.numbers == 1):
	print "[validation] error : argument must be 0 or 1"
	sys.exit(-1)

if  (options.dim > 3 or options.dim < 1):
	print "[validation] error : dimension must be 1, 2 or 3"
	sys.exit(-1)

if not ( options.element == "simplex" or options.element == "hypercube"  ):
	print "[validation] error : element must be simplex or hypercube"
	sys.exit(-1)

if not (options.circle_radius > 0):
	print "[validation] error : radius must be positive"
	sys.exit(-1)

if not (options.line_size > 0):
	print "[validation] error : length must be positive"
	sys.exit(-1)

#print options.order
#print options.numbers
#print options.dim
#print options.element


c = canvas.canvas()


if (options.dim == 1):
        x0 = [0,0]
        x1 = [options.line_size,0]

        hypercube.line(c, x0, x1)

        for i in range(2):
                hypercube.Vertex(c, options.line_size, x0, i, options.circle_radius)
                if ( options.numbers == 1 ):
                        hypercube.numberVertex(c, options.line_size, x0, i, options.circle_radius, 0.1)

        hypercube.Edge(c, options.line_size, x0, 0, options.order, options.circle_radius)

        if ( options.numbers == 1 ):
                hypercube.numberLine(c, options.line_size, x0, options.order, options.circle_radius, 0.1)





if (options.dim == 2):
        x0 = [0,0]

        if ( options.element == "hypercube" ):
                hypercube.border(c, x0, options.line_size)

                for i in range(4):
                        hypercube.Vertex(c, options.line_size, x0, i, options.circle_radius)
                        hypercube.Edge(c, options.line_size, x0, i, options.order, options.circle_radius)

                        if ( options.numbers == 1 ):
                                hypercube.numberVertex(c, options.line_size, x0, i, options.circle_radius, 0.1)
                                hypercube.numberEdge(c, options.line_size, x0, i, options.order, options.circle_radius, 0.1)

                                hypercube.Face(c, options.line_size, x0, options.order, options.circle_radius)

                                if ( options.numbers == 1 ):
                                        hypercube.numberFace(c, options.line_size, x0, options.order, options.circle_radius, 0.1)


        if ( options.element == "simplex" ):
                simplex.border(c, x0, options.line_size)

                for i in range(3):
                        simplex.Vertex(c, options.line_size, x0, i, options.circle_radius)
                        simplex.Edge(c, options.line_size, x0, i, options.order, options.circle_radius)

                        if ( options.numbers == 1 ):
                                simplex.numberVertex(c, options.line_size, x0, i, options.circle_radius, 0.1)
                                simplex.numberEdge(c, options.line_size, x0, i, options.order, options.circle_radius, 0.1)

                        simplex.Face(c, options.line_size, x0, options.order, options.circle_radius)

                        if ( options.numbers == 1 ):
                                simplex.numberFace(c, options.line_size, x0, options.order, options.circle_radius, 0.1)


homename = os.path.expanduser("~/")
dirname = homename + "feel/pointsetvis/equidistributed/" + options.element + "/" + str(options.dim) + "D/P"+ str(options.order)

if not os.path.isdir(dirname + "/"):
        os.makedirs(dirname + "/")

c.writePDFfile(dirname + "/pointset")




# Local Variables:
# indent-tabs-mode: t
# End:

