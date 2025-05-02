import sys


import feelpp.core as feelpp
import feelpp.core.quality as q
import feelpp.toolboxes.core as tb
import feelpp.core.interpolation as I
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.cfpdes import *
import json
import mpi4py
mpi4py.rc.thread_level="single"
import pandas as pd

## Utils 

def remesh_toolbox(f, hclose, hfar, required_facets, required_elts, parent_mesh, cst):
    
        n_required_elts_before=feelpp.nelements(feelpp.markedelements(f.mesh(),required_elts))
        n_required_facets_before=feelpp.nelements(feelpp.markedfaces(f.mesh(),required_facets))
        print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
        print(" . [before remesh] n required facets: {}".format(n_required_facets_before))

        new_mesh, cpt = feelpp.remesh(
            mesh=f.mesh(), metric="gradedls({},{})".format(hclose, hfar), required_elts=required_elts, required_facets=required_facets, params='{"remesh":{ "verbose":-1}}')
        
        print(" . [after remesh]  n remeshes: {}".format(cpt))
        n_required_elts_after=feelpp.nelements(feelpp.markedelements(new_mesh,required_elts))
        n_required_facets_after=feelpp.nelements(feelpp.markedfaces(new_mesh,required_facets))
        print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
        print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
        f.applyRemesh(f.mesh(),new_mesh)

##

## Genral parameters
folder = "shaposwimmer"
sys.argv = [folder]
e = feelpp.Environment(
                sys.argv, opts= feelpp.backend_options("Iv")
                                .add(tb.toolboxes_options("fluid", "fluid"))
                                .add(tb.toolboxes_options("fluid", "pfluid"))
                                .add(tb.toolboxes_options("fluid", "dfluid"))
                                .add(tb.toolboxes_options("coefficient-form-pdes", "expansion")),
                config=feelpp.localRepository(folder)
                                )

hfar = 0.6
hclose = 0.01

required_facets=["Ellipsoid"]
required_elts=["EllipsoidVolume"]

mesh = feelpp.load(feelpp.mesh(dim=3,realdim=3), "test.geo" , 0.1)
feelpp.Environment.setConfigFile('test.cfg')
## Primal problem =====================================================================

#Primal 
fp = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="pfluid", prefix="pfluid")
fp.setMesh(mesh)
fp.init()
fp.printAndSaveInfo()

#remesh_toolbox(f, hclose, hfar, required_facets, required_elts, None, None)

# Reset execution time parameters
fp.reset_executionTime()

#Add Torque
fp.addRigidTorque()
fp.addRigidTorqueRes()

fp.startTimeStep()
     
if feelpp.Environment.isMasterRank():
    print("============================================================\n")
    print("time simulation: {}s iteration : {}\n".format(fp.time(), fp.timeStepBase().iteration()))
    #print("  -- mesh quality: {}s\n".format(min_etaq))
    print("============================================================\n")
    
fp.solve()
fp.exportResults()


## Dual Problem =======================================================================================

## Dual Problem
fd = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="dfluid", prefix="dfluid")
fd.setMesh(mesh)
fd.init()
fd.printAndSaveInfo()


#/!\ faire qu'un seul remaillaige pour avoir un maillage uniforme pour chaque toolboxe
#remesh_toolbox(fd, hclose, hfar, required_facets, required_elts, None, None)


# Reset execution time parameters
fd.reset_executionTime()

#Add Forces
fd.addRigidForce()
fd.addRigidForceRes()

fd.startTimeStep()
     
if feelpp.Environment.isMasterRank():
    print("============================================================\n")
    print("time simulation: {}s iteration : {}\n".format(fd.time(), fd.timeStepBase().iteration()))
    #print("  -- mesh quality: {}s\n".format(min_etaq))
    print("============================================================\n")
    
fd.solve()
fd.exportResults()


## Expansion ===============================================================
exp = cfpdes(dim=3, keyword="expansion", prefix="expansion")
exp.setMesh(mesh)
exp.addParameterInModelProperties("Mu",1.13)
exp.init()
exp.printAndSaveInfo()

exp.solve()
exp.exportResults()
