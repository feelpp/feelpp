import sys


import feelpp.core as feelpp

import feelpp.core.quality as quality
import feelpp.toolboxes.core as tb
import feelpp.core.interpolation as I
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.cfpdes import *
import feelpp.core.meshmover as mm
import mpi4py
mpi4py.rc.thread_level="single"
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

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
                config=feelpp.localRepository(folder))

cst = 0.98
hfar = 0.005
hclose = 0.005
qual = 0.4  #lower bound for the quality of the mesh
al, bl, cl = 0.05, 0.5, 10
ar, br, cr = 0.05, 0.5, 10
l, r = 0, 0



required_facets1=["BoxWalls"]
required_elts1=[]#["EllipsoidVolume"]
required_facets2=["BoxWalls"]
required_elts2=[]#["EllipsoidVolume"]

mesh1 = feelpp.load(feelpp.mesh(dim=3,realdim=3), "fluidandswimmer.geo" , 0.05)
mesh2 = feelpp.load(feelpp.mesh(dim=3,realdim=3), "fluid.geo" , 0.05)


feelpp.Environment.setConfigFile('test2.cfg')
exporter1 = feelpp.exporter(mesh=mesh1, name="fluidandswimmer", geo="change")
exporter2 = feelpp.exporter(mesh=mesh2, name="fluid", geo="change")

List_of_translational_velocity_dot_Text= []
List_of_angular_velocity_dot_Text = []
List_of_remesh_mesh1 = []
List_of_remesh_mesh2 = []
List_of_centermass = []
List_of_volume_swimmer = []

for i in range(200) :
    #interpolation between mesh1 and mesh2
    Pchv1_mesh1 = feelpp.functionSpace(mesh=mesh1, space = "Pchv", order=1)
    Pchv2_mesh1 = feelpp.functionSpace(mesh=mesh1, space = "Pchv", order=2)
    Pchv1_mesh2 = feelpp.functionSpace(mesh=mesh2, space = "Pchv", order=1)
    Pchv2_mesh2 = feelpp.functionSpace(mesh=mesh2, space = "Pchv", order=2)
    interp_laplacian_to_swimmer = I.interpolator(domain = Pchv1_mesh2, image = Pchv1_mesh1,  range = feelpp.elements(mesh1)) 
    interp_swimmer_to_laplacian = I.interpolator(domain = Pchv2_mesh1, image = Pchv2_mesh2,  range = feelpp.elements(mesh2)) 


    ## Dual Problem =======================================================================================

    ## Dual Problem
    fd = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="dfluid", prefix="dfluid")
    fd.setMesh(mesh1)
    fd.init()
    #fd.printAndSaveInfo()


    #/!\ faire qu'un seul remaillaige pour avoir un maillage uniforme pour chaque toolboxe
    #remesh_toolbox(fd, hclose, hfar, required_facets2, required_elts2, None, None)


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
    print("Ud interpolation")
    ud_interp = interp_swimmer_to_laplacian.interpolate(ud)
    print("Ud save")
    fd.saveVelocity(ud_interp, "ud_interp.h5")

    ## Primal problem =====================================================================

    #Primal 
    fp = fluid(dim=3, orderVelocity=2, orderPressure=1, keyword="pfluid", prefix="pfluid")
    fp.setMesh(mesh1)
    fp.init()
    fp.printAndSaveInfo()


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
    up_interp = interp_swimmer_to_laplacian.interpolate(up)
    fp.saveVelocity(up_interp, "up_interp.h5")


    #remesh_toolbox(fp, hclose, hfar, ["Ellipsoid"], ["EllipsoidVolume"], None, None)
    #translationnl velocity
    fp.updateTimeStep()
    print("============================================================\n")
    print("time simulation: ", fp.time(), "s \n")
    print("============================================================\n")
    fp.solve()
    fp.exportResults()

    table_of_fpvalues = pd.read_csv("pfluid.measures/values.csv")
    time_fp = table_of_fpvalues["time"]
    dt_fp = time_fp.values
    dt_fp = dt_fp[0]

    center_of_mass_fp = table_of_fpvalues[["Quantities_body_Ellipsoid.mass_center_0", "Quantities_body_Ellipsoid.mass_center_1", "Quantities_body_Ellipsoid.mass_center_2"]]
    center_of_mass_fp = center_of_mass_fp.values
    translational_velocity_fp = (center_of_mass_fp[0]-center_of_mass_fp[1])/dt_fp
    List_of_translational_velocity_dot_Text.append(np.dot(translational_velocity_fp, np.array([1,0,0])))

    rotation_angle_fp = table_of_fpvalues[["Quantities_body_Ellipsoid.rigid_rotation_angles_0", "Quantities_body_Ellipsoid.rigid_rotation_angles_1", "Quantities_body_Ellipsoid.rigid_rotation_angles_2"]]
    rotation_angle_fp = rotation_angle_fp.values
    angular_velocity_fp = (rotation_angle_fp[0]-rotation_angle_fp[1])/dt_fp
    List_of_angular_velocity_dot_Text.append(np.dot(angular_velocity_fp, np.array([1,0,0])))


 

  

    ## Expansion ===============================================================
    mes = fd.postProcessMeasures().values()
    dfd_mes = pd.DataFrame(mes, index=[0])
    volume_swimmer = dfd_mes["Statistics_volumeswimmer_integrate"][0]
    centermass1 = dfd_mes["Statistics_centermass_integrate_0"][0]/ volume_swimmer
    centermass2 = dfd_mes["Statistics_centermass_integrate_1"][0]/ volume_swimmer
    centermass3 = dfd_mes["Statistics_centermass_integrate_2"][0]/ volume_swimmer
    centermass = np.array([centermass1, centermass2, centermass3])
    Textcrossw = np.cross(np.array([1,0,0]), angular_velocity_fp)

    List_of_volume_swimmer.append(volume_swimmer)
    List_of_centermass.append(centermass)

    
    

    #print("Starting the expansion toolbox")
    exp = cfpdes(dim=3, keyword="expansion", prefix="expansion")
    #print("Setting the mesh for the expansion toolbox")
    exp.setMesh(mesh2)
    print("Initializing the expansion toolbox")
    exp.init()
    print("Adding parameters in the model properties")
    exp.addParameterInModelProperties("Mu",1.13)
    exp.addParameterInModelProperties("volumeswimmer",volume_swimmer)
    exp.addParameterInModelProperties("xCM1", centermass[0])
    exp.addParameterInModelProperties("xCM2", centermass[1])
    exp.addParameterInModelProperties("xCM3", centermass[2])
    exp.addParameterInModelProperties("x01", 0.7)
    exp.addParameterInModelProperties("x02", 0.7)
    exp.addParameterInModelProperties("x03", 0.7)
    exp.addParameterInModelProperties("Textcrossw1", Textcrossw[0])
    exp.addParameterInModelProperties("Textcrossw2", Textcrossw[1])
    exp.addParameterInModelProperties("Textcrossw3", Textcrossw[2])
    exp.addParameterInModelProperties("t",0.005)
    exp.addParameterInModelProperties("l",l)
    exp.addParameterInModelProperties("r", r)
    exp.updateParameterValues()
    #exp.printAndSaveInfo()
    exp.solve()
    exp.exportResults()
    theta = exp.pde("Expansion").fieldUnknown()
    theta_interp = interp_laplacian_to_swimmer.interpolate(theta) 

    mes_exp = exp.postProcessMeasures().values()
    dfexp_mes = pd.DataFrame(mes_exp, index=[0])
    int_grad = dfexp_mes["Statistics_grad_integrate"][0]
    int_surf = dfexp_mes["Statistics_surf_integrate"][0]

    volume0 = np.pi * 0.08**3
    l = al* l + bl*int_grad/int_surf + cl * (volume_swimmer-volume0)/volume0
    r = ar* r + bl*int_grad/int_surf + cr * np.linalg.norm(centermass-np.array([0.7, 0.7, 0.7]))/np.linalg.norm(np.array([0.7, 0.7, 0.7]))

    ## Exporter on the 2 meshes ==========================================================
    exporter1.step(i).setMesh(mesh1)
    exporter1.step(i).add("theta_interp", theta)#_interp)
    exporter1.step(i).add("up", up)
    exporter1.step(i).add("ud", ud)
    exporter1.save()


    exporter2.step(i).setMesh(mesh2)
    exporter2.step(i).add("theta", theta)
    exporter2.step(i).add("up", up_interp)
    exporter2.step(i).add("ud", ud_interp)
    exporter2.save()


    ## Move the meshes ===================================================================
    mesh1 = mm.meshMove(mesh1,theta_interp)
    mesh2 = mm.meshMove(mesh2,theta)

    ## Remesh the meshes ==========================================================
    q1 = quality.etaQ(mesh1).min()
    print(f"q1={q1}")
    if q1 < qual : 
        mesh1, cpt1 = feelpp.remesh(mesh=mesh1, metric=f"gradedls({hclose}, {hfar}, {cst})",required_elts = required_elts1, required_facets=required_facets1,parent=None)
        List_of_remesh_mesh1.append(i)
    print(f"q1 = {quality.etaQ(mesh1).min()}")
    
    q2 = quality.etaQ(mesh2).min()
    print(f"q2={q2}")
    if q2 < qual : 
        mesh2, cpt2 = feelpp.remesh(mesh=mesh2, metric=f"gradedls({hclose}, {hfar}, {cst})",required_elts = required_elts2, required_facets=required_facets2,parent=None)
        List_of_remesh_mesh2.append(i)
    print(f"q2 = {quality.etaQ(mesh2).min()}")



    fig, axs = plt.subplots(3, 2, figsize=(10, 15))
    ax = axs[:,0]


    ax[0].plot(np.array(List_of_translational_velocity_dot_Text)/np.array(List_of_angular_velocity_dot_Text), label=r"$\frac{U\cdot T_{ext}}{\omega\cdot T_{ext}}$")
    ax[1].plot(List_of_translational_velocity_dot_Text, label=r"$U\cdot T_{ext}$")
    ax[2].plot(List_of_angular_velocity_dot_Text, label=r"$\omega\cdot T_{ext}$")

    ax[0].scatter(List_of_remesh_mesh1, np.array(List_of_translational_velocity_dot_Text)[List_of_remesh_mesh1]/np.array(List_of_angular_velocity_dot_Text)[List_of_remesh_mesh1], color='red', marker='o', label='Remesh Mesh1', alpha = 0.8)
    ax[0].scatter(List_of_remesh_mesh2, np.array(List_of_translational_velocity_dot_Text)[List_of_remesh_mesh2]/np.array(List_of_angular_velocity_dot_Text)[List_of_remesh_mesh2], color='green', marker='x', label='Remesh Mesh2')
    ax[1].scatter(List_of_remesh_mesh1, np.array(List_of_translational_velocity_dot_Text)[List_of_remesh_mesh1], color='red', marker='o', label='Remesh Mesh1', alpha = 0.8)
    ax[1].scatter(List_of_remesh_mesh2, np.array(List_of_translational_velocity_dot_Text)[List_of_remesh_mesh2], color='green', marker='x', label='Remesh Mesh2')
    ax[2].scatter(List_of_remesh_mesh1, np.array(List_of_angular_velocity_dot_Text)[List_of_remesh_mesh1], color='red', marker='o', label='Remesh Mesh1', alpha = 0.8)
    ax[2].scatter(List_of_remesh_mesh2, np.array(List_of_angular_velocity_dot_Text)[List_of_remesh_mesh2], color='green', marker='x', label='Remesh Mesh2')
    
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()


    axs[0, 1].plot(np.array(List_of_volume_swimmer)-volume0, label="Error in volume swimmer")
    axs[0, 1].legend()
    axs[0, 1].grid()
    axs[1, 1].plot(np.array(List_of_centermass)[:, 0], label="Centermass x")
    axs[1, 1].plot(np.array(List_of_centermass)[:, 1], label="Centermass y")
    axs[1, 1].plot(np.array(List_of_centermass)[:, 2], label="Centermass z")
    axs[1, 1].legend()
    axs[1, 1].grid()
    axs[2, 1].plot(np.linalg.norm(np.array(List_of_centermass)-np.array([0.5, 0.5, 0.5]), axis=1), label="Error in centermass")
    axs[2, 1].legend()
    axs[2, 1].grid()

    plt.savefig(f"plot.png")
    plt.close(fig)

exporter1.step(i+1).setMesh(mesh1)
exporter1.step(i+1).add("theta_interp", theta_interp)
exporter1.step(i+1).add("up", up)
exporter1.step(i+1).add("ud", ud)
exporter1.save()

exporter2.step(i+1).setMesh(mesh2)
exporter2.step(i+1).add("theta", theta)
exporter2.step(i+1).add("up", up_interp)
exporter2.step(i+1).add("ud", ud_interp)
exporter2.save()