import sys
import feelpp as feelpp
import feelpp.quality as q
import feelpp.toolboxes.core as tb
import feelpp.interpolation as I
from feelpp.toolboxes.fluid import *


sys.argv = ['fluid-remesh']
e = feelpp.Environment(
    sys.argv, opts=tb.toolboxes_options("fluid"),
    config=feelpp.globalRepository("fluid-remesh"))

#    feelpp.Environment.setConfigFile('fluid/TurekHron/cfd3.cfg')
#feelpp.Environment.setConfigFile('fluid/swimmers/3-sphere/2d/three_sphere_2D.cfg')
#feelpp.Environment.setConfigFile('fluid/moving_body/gravity/cfd.cfg')
feelpp.Environment.setConfigFile(
    'cases/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
f = fluid(dim=2, orderVelocity=2, orderPressure=1)
f.init()

hfar = 0.1
hclose = 0.02
def remesh_toolbox(f, hclose, hfar):
    Xh = feelpp.functionSpace(mesh=f.mesh())
    metric = feelpp.gradedls(Xh, feelpp.markedfaces(Xh.mesh(), ["CylinderSurface"]), hclose, hfar)
    R = feelpp.remesher(mesh=f.mesh(), required_elts="CylinderVolume", required_facets="CylinderSurface")
    R.setMetric(metric)
    new_mesh = R.execute()
    f.applyRemesh(new_mesh)

remesh_toolbox(f, hclose, hfar )
f.startTimeStep()
while not f.timeStepBase().isFinished():
    if q.etaQ(f.mesh()).min() < 0.6:
        remesh_toolbox( f, hclose, hfar )
    if feelpp.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: ", f.time(), "s \n")
        print("  -- mesh quality: {}s\n".format(q.etaQ(f.mesh()).min()))
        print("============================================================\n")
    f.solve()
    f.exportResults()
    f.updateTimeStep()
