import feelpp.core as core
import sys

import feelpp.toolboxes.core as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("thermoelectric"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.thermoelectric import *

f=thermoelectric(dim=2,orderDisp=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
