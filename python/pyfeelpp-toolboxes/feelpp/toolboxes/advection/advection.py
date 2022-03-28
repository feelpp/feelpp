from pyfeelpp import core
import sys

import pyfeelpp.toolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("advection"))

#from pyfeelpp import discr,ts,filters
from pyfeelpp.toolboxes.advection import *

f=advection(dim=2,order=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
