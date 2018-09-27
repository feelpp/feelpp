from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("solid"))

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.solid import *

f=solid(dim=2,orderDisp=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
if f.isStationary():
    f.solve()
    f.exportResults()
else:
    if not f.doRestart():
        f.exportResults( f.timeInitial() )
    while not f.timeStepBase().isFinished():
        if f.worldComm().isMasterRank():
            print("============================================================\n")
            print("time simulation: ", f.time(), "s \n")
            print("============================================================\n")
        f.solve()
        f.exportResults()
        f.updateTimeStep()
