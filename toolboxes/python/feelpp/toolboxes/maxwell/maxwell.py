import feelpp.core as core
import sys

import feelpp.toolboxes.core as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("maxwell"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.maxwell import *

f=maxwell(dim=2,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
