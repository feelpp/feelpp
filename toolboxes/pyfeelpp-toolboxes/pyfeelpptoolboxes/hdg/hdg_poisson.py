from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

from pyfeelpptoolboxes.hdg import hdg_mixed_poisson_options
e=core.Environment(sys.argv,opts=hdg_mixed_poisson_options("hdg.mixedpoisson"))

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.hdg import *

f=hdg(dim=2,order=1,worldComm=e.worldCommPtr())
f.init(0,0)
f.printAndSaveInfo()
f.solve()
f.exportResults()
