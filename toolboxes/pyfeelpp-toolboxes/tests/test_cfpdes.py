import sys
import feelpp
import feelpp.toolboxes as tb
import feelpp.toolboxes.cfpdes as cfpdes
import pandas as pd

def test_cfpdes_cfd():
    feelpp.Environment.setConfigFile('cfpdes/fluid/TurekHron/cfd2.cfg')
    f = cfpdes.cfpdes(dim=2)
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
    return not f.checkResults()
#