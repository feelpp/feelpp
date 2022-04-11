import sys
import feelpp
from feelpp.toolboxes.fsi import *
e = feelpp.Environment(sys.argv, opts=toolboxes_options("fsi"))



f=thermoelectric(dim=2,orderDisp=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
