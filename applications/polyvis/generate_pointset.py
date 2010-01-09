from pyx import *
from scipy import *
import math
import square
import sys

from optparse import OptionParser
# Define options
parser = OptionParser()
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

print options.order
print options.numbers

c = canvas.canvas()

x0 = [0,0]

#border of square
line_size = 5
square.border(c, x0, line_size)


for i in range(4):
    square.Vertex(c, line_size, x0, i, 0.07)
    square.Edge(c, line_size, x0, i, options.order, 0.07)

    if ( options.numbers == 1 ):
        print "cucu"
        square.numberVertex(c, line_size, x0, i, 0.07, 0.1)
        square.numberEdge(c, line_size, x0, i, options.order, 0.07, 0.1)

square.Face(c, line_size, x0, options.order, 0.07)

if ( options.numbers == 1 ):
    square.numberFace(c, line_size, x0, options.order, 0.07, 0.1)

c.writePDFfile("pointset_square_" + str(options.order) )
