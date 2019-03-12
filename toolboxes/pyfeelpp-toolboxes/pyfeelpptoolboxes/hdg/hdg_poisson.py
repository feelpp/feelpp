import sys
from pyfeelpptoolboxes.hdg import *

e=core.Environment(sys.argv,opts=hdg_poisson_options())


f=hdgpoisson(dim=2,order=1)
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
