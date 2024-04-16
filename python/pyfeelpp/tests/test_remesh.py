import feelpp.core as fppc
import sys
import pytest
from pathlib import Path
import json

def run(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy = geo
    mesh = fppc.load(m, mesh_name, 0.2)
    
    e = fppc.exporter(mesh=mesh, name=Path(mesh_name).stem, geo="change")
    e.step(0.).setMesh(mesh)
    Xh = fppc.functionSpace(mesh=mesh)
    metric = Xh.element()
    metric.on(range=fppc.elements(mesh), expr=fppc.expr("h/4:h"))
    e.step(0.).add("metric", metric)
    e.save()

    new_mesh, cpt = fppc.remesh(
        mesh=mesh, metric="gradedls({},{})".format(0.02, 0.1), required_elts=[], required_facets=["Gamma_1"], params='{"remesh":{ "verbose":-1}}')

    Xh = fppc.functionSpace(mesh=new_mesh)
    metric = Xh.element()
    metric.on(range=fppc.elements(new_mesh), expr=fppc.expr("h:h"))
    e.step(1.).setMesh(new_mesh)
    e.step(1.).add("metric",metric)
    e.save()
#    assert(new_mesh.numGlobalElements() > mesh.numGlobalElements())


#@pytest.mark.mpi
def test_remesh_2D(init_feelpp):
    
    geo = {
        '2': fppc.create_rectangle()
    }
    # 2D
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/remesh/test_2d")
    run(fppc.mesh(dim=2), geo['2'])

def test_remesh_3D(init_feelpp):
    geo = {
        '3': fppc.create_box()
    }
    # 3D
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/remesh/test_3d")
    run(fppc.mesh(dim=3, realdim=3), geo['3'])
