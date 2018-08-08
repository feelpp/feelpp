from pyfeelpp import *
import sys,time

import pyfeelpp.toolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("solid"))

from pyfeelpp import discr,ts,filters
from pyfeelpp.toolboxes.solid import *

f=Solid_2DP1.create("solid",buildmesh=True,worldComm=e.worldComm())
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
