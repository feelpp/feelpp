import sys

import feelpp.core as fppc
import feelpp.toolboxes as fppt

e=fppc.Environment(sys.argv,opts=tb.toolboxes_options("electric"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.electric import *

f=electric(dim=2,orderPotential=1,worldComm=e.worldCommPtr())
f.init()
#f.printAndSaveInfo()
f.solve()
f.exportResults()
