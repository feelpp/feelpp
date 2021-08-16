from feelpp.toolboxes.fluid import *
import sys
import pytest
import feelpp as feelpp 
import feelpp.toolboxes.core as tb
from feelpp.toolboxes.fluid import *


def test_fluid():
    feelpp.Environment.setConfigFile('fluid/TurekHron/cfd1.cfg')
    # 2D fluid solver using P2P1G1 approximation
    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f.init()
    #f.printAndSaveInfo()
    if f.isStationary():
        f.solve()
        f.exportResults()
    else:
        if not f.doRestart():
            f.exportResults(f.timeInitial())
        while not f.timeStepBase().isFinished():
            if f.worldComm().isMasterRank():
                print("============================================================\n")
                print("time simulation: ", f.time(), "s \n")
                print("============================================================\n")
            f.solve()
            f.exportResults()
            f.updateTimeStep()
