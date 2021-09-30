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
#feelpp.Environment.setConfigFile('cases/moving_body/gravity/cfd.cfg')
#feelpp.Environment.setConfigFile('cases/swimmers/spermatozoon/2d/spermatozoon_2d.cfg')
#feelpp.Environment.setConfigFile('cases/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
f = fluid(dim=2, orderVelocity=2, orderPressure=1)
f.init()

#hfar = 3
#hclose = hfar/10
# 
hfar=0.2
hclose=0.02
# gravity/cfd
hfar=0.2
hclose=hfar/4
# swimmer spermatozzon
hfar=1
hclose=0.75
# 3 spheres
hfar = 3
hclose = hfar/5
def remesh_toolbox(f, hclose, hfar, parent_mesh):
    required_facets=["CircleLeft","CircleCenter","CircleRight"]
    required_elts=["CirLeft","CirCenter","CirRight"]
#    required_facets=["CylinderSurface"]
#    required_elts = ["CylinderVolume"]
    #gravity
    #required_facets=["wall_body"]
    #required_elts = ["Omega_body"]
    #spermatozoon
#    required_facets=["Tail","Head"]
#    required_elts = ["Swimmer"]
    Xh = feelpp.functionSpace(mesh=f.mesh())
    n_required_elts_before=feelpp.nelements(feelpp.markedelements(f.mesh(),required_elts))
    n_required_facets_before=feelpp.nfaces(feelpp.markedfaces(f.mesh(),required_facets))
    print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
    print(" . [before remesh] n required facets: {}".format(n_required_facets_before))
#    metric=feelpp.
    metric = feelpp.gradedls(Xh, feelpp.markedfaces(Xh.mesh(), required_facets), hclose, hfar)
    #R = feelpp.remesher(mesh=f.mesh(), required_elts=required_elts, required_facets=required_facets, parent=parent_mesh)
    R = feelpp.remesher(mesh=f.mesh(), required_elts=required_elts, required_facets=required_facets)
    R.setMetric(metric)
    new_mesh = R.execute()
    n_required_elts_after=feelpp.nelements(feelpp.markedelements(new_mesh,required_elts))
    n_required_facets_after=feelpp.nfaces(feelpp.markedfaces(new_mesh,required_facets))
    print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
    print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
    assert(n_required_elts_before==n_required_elts_after)
    assert(n_required_facets_before==n_required_facets_after)
    f.applyRemesh(new_mesh)
    #assert(feelpp.nfaces(feelpp.markedfaces(f.mesh(), "wall_top"))>0)
    #assert(feelpp.nfaces(feelpp.markedfaces(f.mesh(), "wall_bottom"))>0)
    #assert(feelpp.nfaces(feelpp.markedfaces(f.mesh(), "wall_left"))>0)
    #assert(feelpp.nfaces(feelpp.markedfaces(f.mesh(), "wall_right"))>0)


parent_mesh=f.mesh()
remesh_toolbox(f, hclose, hfar, None )
#remesh_toolbox(f, hclose, hfar, None)
#f.exportResults()
f.startTimeStep()
while not f.timeStepBase().isFinished():
    if q.etaQ(f.mesh()).min() < 0.4:
#    if f.timeStepBase().iteration() % 10 == 0:
        remesh_toolbox(f, hclose, hfar, None)
    if feelpp.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s iteration : {}\n".format(f.time(), f.timeStepBase().iteration()))
        print("  -- mesh quality: {}s\n".format(q.etaQ(f.mesh()).min()))
        print("============================================================\n")
    f.solve()
    f.exportResults()
    
    f.updateTimeStep()
