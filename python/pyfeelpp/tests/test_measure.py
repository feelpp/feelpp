import math
import sys
import feelpp.core as fppc
from fppc.measure import measure
import pytest


def run(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    mesh= fppc.load(m, mesh_name, 0.1)

    M=measure(range=fppc.elements(mesh))
    assert(abs(M-e_meas)<1e-10)
    S_1=measure(range=fppc.markedfaces(mesh,"Gamma_1"))
    assert(abs(S_1-e_s_1) <1e-10)
    S_2 = measure(range=fppc.markedfaces(mesh, "Gamma_2"))
    assert(abs(S_2-e_s_2)<1e-10)
    S_bdy = measure(range=fppc.boundaryfaces(mesh))
    assert(abs(S_bdy-e_s_bdy) < 1e-10)

def test_measure(init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/measure")
    geo = {
        '2': fppc.create_rectangle(),
        '3': fppc.create_box(),
    }
    run(fppc.mesh(dim=2), geo['2'])
    run(fppc.mesh(dim=3, realdim=3), geo['3'])
