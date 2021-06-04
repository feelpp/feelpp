import sys
import feelpp
from feelpp.toolboxes.heat import *

def test_heat():
    feelpp.Environment.setConfigFile(
        'heat/Building/ThermalBridgesENISO10211/thermo2dCase2.cfg')
    f = heat(dim=2, order=1)
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
