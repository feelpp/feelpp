from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("electric"))

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.electric import *

f=electric(dim=2,orderPotential=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
