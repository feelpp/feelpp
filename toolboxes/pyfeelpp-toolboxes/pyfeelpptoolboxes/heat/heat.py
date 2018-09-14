from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("heat"))

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.heat import *

f=heat(dim=2,order=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
