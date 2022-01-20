import sys
import feelpp
from feelpp.toolboxes.core import toolboxes_options

e=feelpp.Environment(sys.argv, opts=toolboxes_options("mixedpoisson", "hdg.poisson"))

from feelpp.toolboxes.hdg import *


f=mixedpoisson(dim=2,order=1)
f.init()
# f.printAndSaveInfo()
f.solve()
f.exportResults()
f.checkResults()
