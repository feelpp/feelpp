import sys
import feelpp.core as fppc
from feelpp.toolboxes.heat import *

e=fppc.Environment(sys.argv,opts=toolboxes_options("heat"))


f=heat(dim=2,order=1)
f.init()
#f.printAndSaveInfo()
f.solve()
f.exportResults()
