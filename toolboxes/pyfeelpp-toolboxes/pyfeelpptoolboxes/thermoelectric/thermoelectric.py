import sys
from pyfeelpptoolboxes.thermoelectric import *

e=core.Environment(sys.argv,opts=toolboxes_options("thermoelectric"))


f=thermoelectric(dim=2,orderDisp=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
