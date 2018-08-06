from pyfeelpp import *
import sys,time

import pyfeelpp.toolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("fluid"))

# from pyfeelpp import discr,ts,filters
from pyfeelpp.toolboxes.fluid import *

f=Fluid_P2P1G1("fluid",buildmesh=True,worldComm=e.worldComm())
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
