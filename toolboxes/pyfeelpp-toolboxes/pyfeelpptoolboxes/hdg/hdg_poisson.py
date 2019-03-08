from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

from pyfeelpptoolboxes.hdg import hdg_poisson_options
e=core.Environment(sys.argv,opts=hdg_poisson_options())

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.hdg import *

f=hdgpoisson(dim=2,order=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
