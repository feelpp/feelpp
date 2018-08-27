from pyfeelpp import core
import sys

import pyfeelpp.toolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("thermoelectric"))

#from pyfeelpp import discr,ts,filters
from pyfeelpp.toolboxes.thermoelectric import *

f=thermoelectric(dim=2,orderDisp=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
