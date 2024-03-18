import math
import sys
import feelpp.core as fppc
from fppc.operators import *
import pytest


def run(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    mesh= fppc.load(m, mesh_name, 0.1)
    Xh = fppc.functionSpace(mesh=mesh)
    v=Xh.element()
    v.on(range=fppc.elements(mesh), expr=fppc.expr("1"))

    M=mass(test=Xh,trial=Xh,range=fppc.elements(mesh))
    assert(abs(M.energy(v,v)-e_meas)<1e-10)
    S=stiffness(test=Xh,trial=Xh,range=fppc.elements(mesh))
    assert(abs(S.energy(v,v)-0)<1e-10)
    v.on(range=fppc.elements(mesh), expr=fppc.expr("x+y:x:y"))
    assert(abs(S.energy(v,v)-2*e_meas)<1e-10)


def test_operators(init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/measure")
    geo = {
        '2': fppc.create_rectangle(),
        '3': fppc.create_box(),
    }
    run(fppc.mesh(dim=2), geo['2'])
    run(fppc.mesh(dim=3, realdim=3), geo['3'])
