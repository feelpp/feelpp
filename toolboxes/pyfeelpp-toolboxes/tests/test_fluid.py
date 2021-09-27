
import sys
import pytest
import feelpp as feelpp 
import feelpp.quality as q
import feelpp.toolboxes.core as tb
import feelpp.interpolation as I
from feelpp.toolboxes.fluid import *


def test_fluid():
    feelpp.Environment.setConfigFile('fluid/TurekHron/cfd3.cfg')
    # 2D fluid solver using P2P1G1 approximation
    f = fluid(dim=2, orderVelocity=2, orderPressure=1,worldComm=feelpp.Environment.worldCommPtr())
    f.init()
    #f.printAndSaveInfo()
    if f.isStationary():
        f.solve()
        f.exportResults()
    else:
        if not f.doRestart():
            f.exportResults(f.timeInitial())
        f.startTimeStep()
        while not f.timeStepBase().isFinished():
            if feelpp.Environment.isMasterRank(): #.worldComm().isMasterRank():
                print("============================================================\n")
                print("time simulation: ", f.time(), " - ", "s \n")
                print("============================================================\n")
            f.solve()
            f.exportResults()
            f.updateTimeStep()


def test_fluid_remesh():
    feelpp.Environment.setConfigFile('fluid/TurekHron/cfd3.cfg')
    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f.init()

    e = feelpp.exporter(mesh=f.mesh(), name="turek-hron", geo="change")
    e.step(0.).setMesh(f.mesh())
    e.step(0.).add("velocity",f.fieldVelocity())
    e.step(0.).add("pressure", f.fieldPressure())
    e.save()
    f.startTimeStep()
    while f.time() < f.timeFinal()/2:
        if feelpp.Environment.isMasterRank():
            print("============================================================\n")
            print("time simulation: ", f.time(), "s \n")
            print("============================================================\n")

        f.solve()
        e.step(f.time()).add("velocity", f.fieldVelocity())
        e.step(f.time()).add("pressure", f.fieldPressure())
        e.save()
        f.updateTimeStep()

    hclose=0.005
    hfar=0.07
    Xh = feelpp.functionSpace(mesh=f.mesh())
    metric = feelpp.gradedls(Xh, feelpp.markedfaces(Xh.mesh(),["wall2","wake"]), hclose, hfar)

    R = feelpp.remesher(mesh=f.mesh())
    R.setMetric(metric)
    new_mesh = R.execute()

    f2 = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f2.setMesh(new_mesh)
    f2.setTimeInitial(f.time()-f2.timeStep())
    f2.init()
    f2.init(f,["Fluid"])
    f2.setTimeInitial(f.time()-f2.timeStep())
    #OIv = I.interpolator(domain=f.functionSpaceVelocity(), image=f2.functionSpaceVelocity(), range=feelpp.elements(f2.mesh()))
    #f2.setFieldVelocity(OIv.interpolate(domain=f.fieldVelocity()))
    #OIp = I.interpolator(domain=f.functionSpacePressure(), image=f2.functionSpacePressure(), range=feelpp.elements(f2.mesh()))
    #f2.setFieldPressure(OIp.interpolate(domain=f.fieldPressure()))
    print("time: {} step: {} initial: {}".format(f2.time(), f2.timeStep(), f2.timeInitial()))
    
    #e.step(f2.time()+f2.timeStep()).setMesh(f2.mesh())
    #e.step(f2.time()+f2.timeStep()).add("velocity", f2.fieldVelocity())
    #e.step(f2.time()+f2.timeStep()).add("pressure", f2.fieldPressure())
    #e.save()
    #return
    f2.startTimeStep()

    print("time: {} step: {} initial: {}".format(
        f2.time(), f2.timeStep(), f2.timeInitial()))
    while not f2.timeStepBase().isFinished():
        if feelpp.Environment.isMasterRank():
            print("============================================================\n")
            print("time simulation: ", f2.time(), "s \n")
            print("============================================================\n")
        f2.solve()
        e.step(f2.time()).setMesh(f2.mesh())
        e.step(f2.time()).add("velocity",f2.fieldVelocity())
        e.step(f2.time()).add("pressure", f2.fieldPressure())
        e.save()
        f2.updateTimeStep()
    
