import sys

import feelpp.toolboxes.modelcore as modelcore

e=feelpp.Environment(sys.argv,opts=modelcore.toolboxes_options("electric"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.electric import *

f=electric(dim=2,orderPotential=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
