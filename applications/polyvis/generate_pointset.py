from pyx import *
from scipy import *

import os
import sys
import string

import math
import square
import triangle


from optparse import OptionParser
# Define options
parser = OptionParser()
parser.add_option("--dim", dest="dim", type="int", help="topological dimension of the element", default=2)
parser.add_option("--element_type", dest="element", type="string", help="type of the reference element: simplex or hypercube", default="simplex")
parser.add_option("--order", dest="order", type="int", help="pointset order", default=1)
parser.add_option("--numbering", dest="numbers", type="int", help="adds numbering to points: 0 for no numbering, 1 to add the numbers", default=0)
(options, args) = parser.parse_args()


# Verification of the options
if (options.order < 0):
	print "[validation] error : order must be nonnegative"
	sys.exit(-1)

if (options.numbers < 0):
	print "[validation] error : argument must be 0 or 1"
	sys.exit(-1)

if (options.numbers > 1):
	print "[validation] error : argument must be 0 or 1"
	sys.exit(-1)

if (options.dim < 1):
	print "[validation] error : argument must be 1, 2 or 3"
	sys.exit(-1)

if (options.dim > 3):
	print "[validation] error : argument must be 1, 2 or 3"
	sys.exit(-1)

print options.order
print options.numbers
print options.dim
print options.element


c = canvas.canvas()


if (options.dim == 2):
        x0 = [0,0]

        line_size = 5

        if ( options.element == "hypercube" ):
                square.border(c, x0, line_size)

                for i in range(4):
                        square.Vertex(c, line_size, x0, i, 0.07)
                        square.Edge(c, line_size, x0, i, options.order, 0.07)

                        if ( options.numbers == 1 ):
                                square.numberVertex(c, line_size, x0, i, 0.07, 0.1)
                                square.numberEdge(c, line_size, x0, i, options.order, 0.07, 0.1)

                                square.Face(c, line_size, x0, options.order, 0.07)

                                if ( options.numbers == 1 ):
                                        square.numberFace(c, line_size, x0, options.order, 0.07, 0.1)


        if ( options.element == "simplex" ):
                triangle.border(c, x0, line_size)

                for i in range(3):
                        triangle.Vertex(c, line_size, x0, i, 0.07)
                        triangle.Edge(c, line_size, x0, i, options.order, 0.07)

                        if ( options.numbers == 1 ):
                                triangle.numberVertex(c, line_size, x0, i, 0.07, 0.1)
                                triangle.numberEdge(c, line_size, x0, i, options.order, 0.07, 0.1)

                        triangle.Face(c, line_size, x0, options.order, 0.07)

                        if ( options.numbers == 1 ):
                                triangle.numberFace(c, line_size, x0, options.order, 0.07, 0.1)


homename = os.path.expanduser("~/")
dirname = homename + "life/pointsetvis/equidistributed/" + options.element + "/" + str(options.dim) + "D/P"+ str(options.order)

print dirname

if not os.path.isdir(dirname + "/"):
        os.makedirs(dirname + "/")

c.writePDFfile(dirname + "/pointset")
