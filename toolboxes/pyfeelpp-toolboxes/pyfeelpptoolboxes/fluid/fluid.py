from pyfeelpp import core
import sys

import pyfeelpptoolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("fluid"))

# from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.fluid import *

# 2D fluid solver using P2P1G1 approximation
f=fluid(dim=2,orderVelocity=2,orderPressure=1,
        worldComm=core.Environment.worldCommPtr())
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
