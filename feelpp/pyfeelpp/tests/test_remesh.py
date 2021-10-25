import feelpp
import sys
import pytest
from pathlib import Path

def run(m, geo):
    mesh_name, e_meas, e_s_1, e_s_2, e_s_bdy = geo
    mesh = feelpp.load(m, mesh_name, 0.2)
    
    e = feelpp.exporter(mesh=mesh, name=Path(mesh_name).stem, geo="change")
    e.step(0.).setMesh(mesh)
    Xh = feelpp.functionSpace(mesh=mesh)
    metric = Xh.element()
    metric.on(range=feelpp.elements(mesh), expr=feelpp.expr("h/4:h"))
    e.step(0.).add("metric", metric)
    e.save()

    new_mesh, cpt = feelpp.remesh(
        mesh, "gradedls({},{})".format(0.02, 0.1), [], ["Gamma_1"], None)

    Xh = feelpp.functionSpace(mesh=new_mesh)
    metric = Xh.element()
    metric.on(range=feelpp.elements(new_mesh), expr=feelpp.expr("h:h"))
    e.step(1.).setMesh(new_mesh)
    e.step(1.).add("metric",metric)
    e.save()
#    assert(new_mesh.numGlobalElements() > mesh.numGlobalElements())


#@pytest.mark.mpi
def test_remesh_2D(init_feelpp):
    
    geo = {
        '2': feelpp.create_rectangle()
    }
    # 2D
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/remesh/test_2d")
    run(feelpp.mesh(dim=2), geo['2'])

def test_remesh_3D(init_feelpp):
    geo = {
        '3': feelpp.create_box()
    }
    # 3D
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/remesh/test_3d")
    run(feelpp.mesh(dim=3, realdim=3), geo['3'])
