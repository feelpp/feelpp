import sys
import feelpp
from feelpp.toolboxes.heat import *

e=feelpp.Environment(sys.argv,opts=toolboxes_options("heat"))


f=heat(dim=2,order=1)
f.init()
#f.printAndSaveInfo()
f.solve()
f.exportResults()
