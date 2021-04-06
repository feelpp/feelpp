import sys
from feelpp.toolboxes.thermoelectric import *

e=feelpp.Environment(sys.argv,opts=toolboxes_options("thermo-electric"))


f=thermoelectric(dim=2,orderPotential=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
