
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
#    feelpp.Environment.setConfigFile('fluid/TurekHron/cfd3.cfg')
    #feelpp.Environment.setConfigFile('fluid/swimmers/3-sphere/2d/three_sphere_2D.cfg')
    #feelpp.Environment.setConfigFile('fluid/moving_body/gravity/cfd.cfg')
    feelpp.Environment.setConfigFile(
        'fluid/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
    f.init()

    e = feelpp.exporter(mesh=f.mesh(), name="turek-hron", geo="change")
    e.step(0.).setMesh(f.mesh())
    e.step(0.).add("velocity",f.fieldVelocity())
    e.step(0.).add("pressure", f.fieldPressure())
    e.save()
    f.startTimeStep()
    while not f.timeStepBase().isFinished():

        if f.timeStepBase().iteration() % 4 == 0:
            hfar=0.1
            hclose=0.02
            Xh = feelpp.functionSpace(mesh=f.mesh())
            metric = feelpp.gradedls(Xh, feelpp.markedfaces(
                Xh.mesh(), ["CylinderSurface"]), hclose, hfar)
            R = feelpp.remesher(mesh=f.mesh(),required_elts="CylinderVolume",required_facets="CylinderSurface")
            R.setMetric(metric)
            new_mesh = R.execute()
            f.applyRemesh(new_mesh)

        if feelpp.Environment.isMasterRank():
            print("============================================================\n")
            print("time simulation: ", f.time(), "s \n")
            print("============================================================\n")

        f.solve()
        e.step(f.time()).setMesh(f.mesh())
        e.step(f.time()).add("velocity", f.fieldVelocity())
        e.step(f.time()).add("pressure", f.fieldPressure())
        e.save()
        f.updateTimeStep()

    
