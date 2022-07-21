import math
import sys
import feelpp
from feelpp.integrate  import integrate
import pytest

def run(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy = geo
    mesh = feelpp.load(m, mesh_name, 0.1)

    M = integrate(range=feelpp.elements(mesh),expr="1")
    assert(abs(M[0]-e_meas) < 1e-10)
    S_1 = integrate(range=feelpp.markedfaces(mesh, "Gamma_1"), expr="1")
    assert(abs(S_1[0]-e_s_1) < 1e-10)
    S_2 = integrate(range=feelpp.markedfaces(mesh, "Gamma_2"), expr="1")
    assert(abs(S_2[0]-e_s_2) < 1e-10)
    S_bdy = integrate(range=feelpp.boundaryfaces(mesh), expr="1")
    assert(abs(S_bdy[0]-e_s_bdy) < 1e-10)
    if mesh.dimension() == 2 :
        S_bdy = integrate(range=feelpp.boundaryfaces(mesh), expr="nx+ny:nx:ny")
    elif mesh.dimension() == 3:
        S_bdy = integrate(range=feelpp.boundaryfaces(mesh), expr="nx+ny+nz:nx:ny:nz")
    assert(abs(S_bdy[0]) < 1e-10)

def test_integrate(init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/integrate")
    geo = {
        '2': feelpp.create_rectangle(),
        '3': feelpp.create_box(),
    }
    run(feelpp.mesh(dim=2), geo['2'])
    run(feelpp.mesh(dim=3, realdim=3), geo['3'])
