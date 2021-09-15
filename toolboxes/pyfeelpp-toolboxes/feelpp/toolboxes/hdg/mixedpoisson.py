import sys
import feelpp
from feelpp.toolboxes.hdg import *

e=feelpp.Environment(sys.argv, opts=toolboxes_options("mixedpoisson", "hdg.poisson"))

f=mixedpoisson(dim=2,order=1)
f.init()
# f.printAndSaveInfo()
f.solve()
f.exportResults()

