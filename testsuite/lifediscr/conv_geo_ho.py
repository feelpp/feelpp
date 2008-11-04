#! /usr/bin/python
import os
import sys
import re
from pyx import *
from pyx.deco import barrow,earrow
from pyx.style import linewidth, linestyle
from pyx.graph.axis import painter, tick
from pyx.graph.axis import *
from scipy import *
from scipy.optimize import leastsq
from optparse import OptionParser
from subprocess import Popen, PIPE

parser = OptionParser()
parser.add_option("--ntests", type="int",dest="ntests", help="number of tests", default=3)
parser.add_option("--hsize", type="float",dest="hsize", help="starting h size", default=.5)
(options, args) = parser.parse_args()

v = zeros([options.ntests,6])
print v
for i in range(1,options.ntests+1) :
    for order in range(1,5+1) :
        print 'Execute case h='+str(options.hsize/i)+' Order='+str(order)
        v[i-1][0]=options.hsize/i
        proc = Popen(['test_integration_ho', '--order='+str(order),'--hsize='+str(options.hsize/i)],stdout=PIPE, stderr=PIPE)
        return_code = proc.wait()
        if return_code == 0:
            std=proc.stdout.read()
            print "Failure %s:\n%s" % (return_code, std)

            res=re.search('int\( 1 \)=([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)',std).group(1)
            print "  integral = %s" % res
            v[i-1][order]=abs(pi-float(res))
        else:
            print "Failure %s:\n%s" % (return_code, proc.stderr.read())

io.write_array("conv.dat", v, precision=16)
