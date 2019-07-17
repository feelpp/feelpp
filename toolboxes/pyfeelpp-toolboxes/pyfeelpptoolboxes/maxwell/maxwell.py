from pyfeelpp import core
import sys

import pyfeelpp.toolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("maxwell"))

#from pyfeelpp import discr,ts,filters
from pyfeelpp.toolboxes.maxwell import *

f=maxwell(dim=2,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
