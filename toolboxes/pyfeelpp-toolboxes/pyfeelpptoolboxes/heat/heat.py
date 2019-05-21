import sys
from pyfeelpptoolboxes.heat import *

e=core.Environment(sys.argv,opts=toolboxes_options("heat"))


f=heat(dim=2,order=1)
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
