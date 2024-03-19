import sys
import feelpp.core as fppc
from feelpp.toolboxes.solid import *
e = fppc.Environment(sys.argv, opts=toolboxes_options("solid"))


f=solid(dim=2,orderDisp=1)
f.init()
#f.printAndSaveInfo()
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
