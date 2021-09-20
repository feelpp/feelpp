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
    e.step(0.).add("metric",metric)

    e.save()

    R = feelpp.remesher(mesh=mesh)
    R.setMetric(metric)
    new_mesh=R.execute()
    e.step(1.).setMesh(new_mesh)
    Xh = feelpp.functionSpace(mesh=new_mesh)
    metric = Xh.element()
    metric.on(range=feelpp.elements(mesh), expr=feelpp.expr("h:h"))
    e.step(1.).add("metric",metric)
    e.save()
    assert(new_mesh.numGlobalElements() > mesh.numGlobalElements())


#@pytest.mark.mpi
def test_remesh(init_feelpp):
    
    geo = {
        '2': feelpp.create_rectangle(),
        '3': feelpp.create_box()
    }
    # no parallel remeshing in 2D
    if feelpp.Environment.numberOfProcessors() == 1:
        feelpp.Environment.changeRepository(
            directory="pyfeelpp-tests/remesh/test_2d")
        run(feelpp.mesh(dim=2), geo['2'])
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/remesh/test_3d")
    run(feelpp.mesh(dim=3, realdim=3), geo['3'])
