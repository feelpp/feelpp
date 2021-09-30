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
feelpp.Environment.setConfigFile('cases/swimmers/3-sphere/2d/three_sphere_2D.cfg')
#feelpp.Environment.setConfigFile('fluid/moving_body/gravity/cfd.cfg')
#feelpp.Environment.setConfigFile(
#    'cases/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
f = fluid(dim=2, orderVelocity=2, orderPressure=1)
f.init()

hfar = 3
hclose = hfar/10
def remesh_toolbox(f, hclose, hfar, parent_mesh):
    required_facets=["CircleLeft","CircleCenter","CircleRight"]
    required_elts=["CirLeft","CirCenter","CirRight"]
    Xh = feelpp.functionSpace(mesh=f.mesh())
    n_required_elts_before=feelpp.nelements(feelpp.markedelements(f.mesh(),required_elts))
    n_required_facets_before=feelpp.nfaces(feelpp.markedfaces(f.mesh(),required_facets))
    print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
    print(" . [before remesh] n required facets: {}".format(n_required_facets_before))
    metric = feelpp.gradedls(Xh, feelpp.markedfaces(Xh.mesh(), required_facets), hclose, hfar)
    R = feelpp.remesher(mesh=f.mesh(), required_elts=required_elts, required_facets=required_facets, parent=parent_mesh)
    R.setMetric(metric)
    new_mesh = R.execute()
    n_required_elts_after=feelpp.nelements(feelpp.markedelements(new_mesh,required_elts))
    n_required_facets_after=feelpp.nfaces(feelpp.markedfaces(new_mesh,required_facets))
    print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
    print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
    assert(n_required_elts_before==n_required_elts_after)
    assert(n_required_facets_before==n_required_facets_after)
    f.applyRemesh(new_mesh)

parent_mesh=f.mesh()
remesh_toolbox(f, hclose, hfar )
f.startTimeStep()
while not f.timeStepBase().isFinished():
#    if q.etaQ(f.mesh()).min() < 0.6:
    if f.timeStepBase().iteration()%1 == 0:
        remesh_toolbox( f, hclose, hfar, parent_mesh )
    if feelpp.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: ", f.time(), "s \n")
        print("  -- mesh quality: {}s\n".format(q.etaQ(f.mesh()).min()))
        print("============================================================\n")
    f.solve()
    f.exportResults()
    f.updateTimeStep()
