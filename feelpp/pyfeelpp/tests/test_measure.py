import math
import sys
import feelpp
import pytest


def create_rectangle():
    if feelpp.Environment.isMasterRank():
        import gmsh
        gmsh.initialize()
        gmsh.model.add("rectangle")
        gmsh.logger.start()
        rectangle=gmsh.model.occ.addRectangle(0, 0, 0, 1, 2)
        gmsh.model.occ.synchronize()
        omega=gmsh.model.addPhysicalGroup(2, [rectangle])
        gmsh.model.setPhysicalName(2,omega,"Omega")
        gamma_1=gmsh.model.addPhysicalGroup(1, [1, 3])
        gmsh.model.setPhysicalName(1, gamma_1, "Gamma_1")
        gamma_2=gmsh.model.addPhysicalGroup(1, [2, 4])
        gmsh.model.setPhysicalName(1, gamma_2, "Gamma_2")
        gmsh.model.mesh.generate(2)
        gmsh.write("rectangle.msh")
    return "rectangle.msh",2,2,4,6


def create_box():
    if feelpp.Environment.isMasterRank():
        import gmsh
        gmsh.initialize()
        gmsh.model.add("box")
        gmsh.logger.start()
        box = gmsh.model.occ.addBox(0, 0, 0, 1, 2, 0.5)
        gmsh.model.occ.synchronize()
        omega = gmsh.model.addPhysicalGroup(3, [box])
        gmsh.model.setPhysicalName(3, omega, "Omega")
        gamma_1 = gmsh.model.addPhysicalGroup(2, [1, 3])
        gmsh.model.setPhysicalName(2, gamma_1, "Gamma_1")
        gamma_2 = gmsh.model.addPhysicalGroup(2, [2, 4])
        gmsh.model.setPhysicalName(2, gamma_2, "Gamma_2")
        gamma_3= gmsh.model.addPhysicalGroup(2, [5, 6])
        gmsh.model.setPhysicalName(2, gamma_3, "Gamma_3")
        gmsh.model.mesh.generate(3)
        gmsh.write("box.msh")
    return "box.msh", 1, 1.5, 1.5,7



def run(m, geo):
    mesh_name, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    mesh= feelpp.load(m, mesh_name, 0.1)

    M=feelpp.measure(range=feelpp.elements(mesh))
    assert(abs(M-e_meas)<1e-10)
    S_1=feelpp.measure(range=feelpp.markedfaces(mesh,"Gamma_1"))
    assert(abs(S_1-e_s_1) <1e-10)
    S_2 = feelpp.measure(range=feelpp.markedfaces(mesh, "Gamma_2"))
    assert(abs(S_2-e_s_2)<1e-10)
    S_bdy = feelpp.measure(range=feelpp.boundaryfaces(mesh))
    assert(abs(S_bdy-e_s_bdy) < 1e-10)

def test_measure(init_feelpp):
    geo = {
        '2': create_rectangle(),
        '3': create_box(),
    }
    run(feelpp.mesh(dim=2), geo['2'])
    run(feelpp.mesh(dim=3, realdim=3), geo['3'])
