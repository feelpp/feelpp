import feelpp.core as core
import sys

import feelpp.toolboxes.core as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("advection"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.advection import *

f=advection(dim=2,order=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
