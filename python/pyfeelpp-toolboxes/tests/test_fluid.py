
import sys
import pytest
import feelpp.core as fppc
import feelpp.core.interpolation as I
import feelpp.core.quality as q
import feelpp.toolboxes.core as tb
from feelpp.toolboxes.fluid import *


#@pytest.mark.order("first")

def test_fluid1():
    fppc.Environment.setConfigFile('fluid/TurekHron/cfd1.cfg')
    # 2D fluid solver using P2P1G1 approximation
    f = fluid(dim=2, orderVelocity=2, orderPressure=1,worldComm=fppc.Environment.worldCommPtr())
    simulate(f)

#@pytest.mark.order("second")
#def test_fluid2_remesh():
##    fppc.Environment.setConfigFile('fluid/TurekHron/cfd3.cfg')
#    #fppc.Environment.setConfigFile('fluid/swimmers/3-sphere/2d/three_sphere_2D.cfg')
#    #fppc.Environment.setConfigFile('fluid/moving_body/gravity/cfd.cfg')
#    fppc.Environment.setConfigFile(
#        'fluid/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
#    f = fluid(dim=2, orderVelocity=2, orderPressure=1)
#    f.init()
#
#    f.exportResults()
#    f.startTimeStep()
#    while not f.timeStepBase().isFinished():
#
#        if f.timeStepBase().iteration() % 4 == 0:
#            hfar=0.1
#            hclose=0.02
#            Xh = fppc.functionSpace(mesh=f.mesh())
#            metric = fppc.gradedls(Xh, fppc.markedfaces(
#                Xh.mesh(), ["CylinderSurface"]), hclose, hfar)
#            R = fppc.remesher(mesh=f.mesh(),required_elts="CylinderVolume",required_facets="CylinderSurface")
#            R.setMetric(metric)
#            new_mesh = R.execute()
#            f.applyRemesh(new_mesh)
#
#        if fppc.Environment.isMasterRank():
#            print("============================================================\n")
#            print("time simulation: ", f.time(), "s \n")
#            print("============================================================\n")
#
#        f.solve()
#        f.exportResults()
#        f.updateTimeStep()
#

    
