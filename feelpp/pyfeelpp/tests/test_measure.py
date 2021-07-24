import math
import sys
import feelpp
from feelpp.measure import measure
import pytest


def run(m, geo):
    mesh_name, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    mesh= feelpp.load(m, mesh_name, 0.1)

    M=measure(range=feelpp.elements(mesh))
    assert(abs(M-e_meas)<1e-10)
    S_1=measure(range=feelpp.markedfaces(mesh,"Gamma_1"))
    assert(abs(S_1-e_s_1) <1e-10)
    S_2 = measure(range=feelpp.markedfaces(mesh, "Gamma_2"))
    assert(abs(S_2-e_s_2)<1e-10)
    S_bdy = measure(range=feelpp.boundaryfaces(mesh))
    assert(abs(S_bdy-e_s_bdy) < 1e-10)

def test_measure(init_feelpp):
    geo = {
        '2': feelpp.create_rectangle(),
        '3': feelpp.create_box(),
    }
    run(feelpp.mesh(dim=2), geo['2'])
    run(feelpp.mesh(dim=3, realdim=3), geo['3'])
