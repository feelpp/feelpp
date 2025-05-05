import sys


import feelpp.core as feelpp
import feelpp.core.quality as q
import feelpp.toolboxes.core as tb
import feelpp.core.interpolation as I
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.cfpdes import *
import feelpp.core.meshmover as mm
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

mesh1 = feelpp.load(feelpp.mesh(dim=3,realdim=3), "fluidandswimmer.geo" , 0.1)
mesh2 = feelpp.load(feelpp.mesh(dim=3,realdim=3), "fluid.geo" , 0.1)
#interpolation between mesh1 and mesh2
Pchv1_mesh1 = feelpp.functionSpace(mesh=mesh1, space = "Pchv", order=1)
Pchv1_mesh2 = feelpp.functionSpace(mesh=mesh2, space = "Pchv", order=1)
# interp = I.interpolator(domain = Pchv1_mesh2, image = Pchv1_mesh1,  range = feelpp.elements(mesh1)) PCHV3D PAS PRIS EN COMPTE


feelpp.Environment.setConfigFile('test.cfg')
exporter1 = feelpp.exporter(mesh=mesh1, name="fluidandswimmer", geo="change")
exporter2 = feelpp.exporter(mesh=mesh2, name="fluid", geo="change")
## Primal problem =====================================================================

for i in range(1) :

    #Primal 
    fp = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="pfluid", prefix="pfluid")
    fp.setMesh(mesh1)
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
    up = fp.fieldVelocity()


    ## Dual Problem =======================================================================================

    ## Dual Problem
    fd = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="dfluid", prefix="dfluid")
    fd.setMesh(mesh1)
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
    ud = fd.fieldVelocity()


    ## Expansion ===============================================================
    exp = cfpdes(dim=3, keyword="expansion", prefix="expansion")
    exp.setMesh(mesh2)
    exp.addParameterInModelProperties("Mu",1.13)
    exp.addParameterInModelProperties("t",1.)
    exp.addParameterInModelProperties("l",0)
    exp.init()
    exp.printAndSaveInfo()

    exp.solve()
    exp.exportResults()
    theta = exp.pde("Expansion").fieldUnknown()
    
    print(dir(theta))
    #theta_interp = interp.interpolate(theta) PCHV3D PAS PRIS EN COMPTE

    #thetax = theta.comp("X") NE MARCHE PAS EN PYTHON
    #/!\ changer espace d'interpolation en Pch !
    # theta_interp_x = interp.interpolate(xtheta)
    # theta_interp_y = interp.interpolate(ytheta)
    # theta_interp_z = interp.interpolate(ztheta)

    exporter1.step(i).setMesh(mesh1)
    exporter1.step(i).add("theta", theta)
    exporter1.step(i).add("up", up)
    exporter1.step(i).add("ud", ud)
    exporter1.save()


    exporter2.step(i).setMesh(mesh2)
    exporter2.step(i).add("theta", theta)
    exporter2.step(i).add("up", up)
    exporter2.step(i).add("ud", ud)
    exporter2.save()


    #mesh1 = mm.meshMove(mesh1,theta)
    mesh2 = mm.meshMove(mesh2,theta)


    

exporter1.step(i+1).setMesh(mesh1)
exporter1.step(i+1).add("theta", theta)
exporter1.step(i+1).add("up", up)
exporter1.step(i+1).add("ud", ud)
exporter1.save()

exporter2.step(i+1).setMesh(mesh2)
exporter2.step(i+1).add("theta", theta)
exporter2.step(i+1).add("up", up)
exporter2.step(i+1).add("ud", ud)
exporter2.save()