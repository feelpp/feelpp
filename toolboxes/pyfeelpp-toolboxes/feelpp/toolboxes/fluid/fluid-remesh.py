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
#feelpp.Environment.setConfigFile('cases/swimmers/3-sphere/2d/three_sphere_2D.cfg')
#feelpp.Environment.setConfigFile('cases/swimmers/3-sphere/2d-shear/three_sphere_2D_shearflow.cfg')
#feelpp.Environment.setConfigFile('cases/moving_body/gravity/cfd.cfg')
#feelpp.Environment.setConfigFile('cases/swimmers/spermatozoon/2d/spermatozoon_2d.cfg')
#feelpp.Environment.setConfigFile('cases/moving_body/gravity/cylinder_under_gravity/cylinder_under_gravity.cfg')
feelpp.Environment.setConfigFile('cases/moving_body/gravity/shape_under_gravity/2d/shape_under_gravity.cfg')
#feelpp.Environment.setConfigFile('cases/moving_body/gravity/shape_under_gravity/3d/shape_under_gravity.cfg')
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
## gravity 3D
hfar=1.
#hclose=0.1
# swimmer spermatozzon
#hfar=1
#hclose=0.75
# 3 spheres
#hfar = 1
#hclose = 0.002
def remesh_toolbox(f, hclose, hfar, parent_mesh):
    # 3 spheres
    #required_facets=["CircleLeft","CircleCenter","CircleRight"]
    #required_elts=["CirLeft","CirCenter","CirRight"]
#    required_facets=["CylinderSurface"]
#    required_elts = ["CylinderVolume"]
    #gravity
    #required_facets=["wall_body"]
    #required_elts = ["Omega_body"]
    required_facets=["Wall_Body"]
    required_elts = ["Body"]
    #spermatozoon
#    required_facets=["Tail","Head"]
#    required_elts = ["Swimmer"]
    Xh = feelpp.functionSpace(mesh=f.mesh())
    n_required_elts_before=feelpp.nelements(feelpp.markedelements(f.mesh(),required_elts))
    n_required_facets_before=feelpp.nfaces(feelpp.markedfaces(f.mesh(),required_facets))
    print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
    print(" . [before remesh] n required facets: {}".format(n_required_facets_before))
#    metric=feelpp.
    new_mesh,cpt = feelpp.remesh(f.mesh(), "gradedls({},{})".format(hclose, hfar), required_elts, required_facets, None )
    print(" . [after remesh]  n remeshes: {}".format(cpt))
    n_required_elts_after=feelpp.nelements(feelpp.markedelements(new_mesh,required_elts))
    n_required_facets_after=feelpp.nfaces(feelpp.markedfaces(new_mesh,required_facets))
    print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
    print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
    f.applyRemesh(new_mesh)

parent_mesh=f.mesh()
#remesh_toolbox(f, hclose, hfar, None )
#remesh_toolbox(f, hclose, hfar, None)
#f.exportResults()
f.startTimeStep()
while not f.timeStepBase().isFinished():
    min_etaq = q.etaQ(f.mesh()).min()
    if min_etaq < 0.72:
#    if f.timeStepBase().iteration() % 10 == 0:
        remesh_toolbox(f, hclose, hfar, None)
    if feelpp.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s iteration : {}\n".format(f.time(), f.timeStepBase().iteration()))
        print("  -- mesh quality: {}s\n".format(min_etaq))
        print("============================================================\n")
    for i in range(4):
        f.solve()
    f.exportResults()
    
    f.updateTimeStep()
