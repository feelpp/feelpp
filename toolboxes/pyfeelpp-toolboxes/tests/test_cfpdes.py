import sys
import feelpp
import feelpp.toolboxes.core as tb
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


def test_cfpdes_remesh():
    feelpp.Environment.setConfigFile('cfpdes/laplace/l-shape/l-shape.cfg')
    f = cfpdes.cfpdes(dim=2)
    f.init()
    f.solve()
    e = feelpp.exporter(mesh=f.mesh(), name="l-shape", geo="change")
    e.step(0.).setMesh(f.mesh())
    f.exportSolutionToStep( e.step(0.) )
    e.save()
    Xh = feelpp.functionSpace(mesh=f.mesh())
    metric = Xh.element()
    metric.on(range=feelpp.elements(f.mesh()), expr=feelpp.expr("0.1"))
    R = feelpp.remesher(mesh=f.mesh())
    R.setMetric(metric)
    new_mesh = R.execute()
    
    
    fnew = cfpdes.cfpdes(dim=2)
    fnew.setMesh(new_mesh)
    fnew.init()
    fnew.solve()
    e.step(1.).setMesh(new_mesh)
    fnew.exportSolutionToStep(e.step(1.))
    e.save()
    return not fnew.checkResults()
