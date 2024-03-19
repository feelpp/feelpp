import feelpp.core as fppc
import sys,os
import pytest
from feelpp.core.timing import tic,toc

cases = [
         (2,'feelpp2d',os.path.dirname(__file__)+'/cases/feelpp2d/feelpp2d.geo'),
         (3,'feelpp3d',os.path.dirname(__file__)+'/cases/feelpp3d/feelpp3d.geo')
        ]


#@pytest.mark.mpi
@pytest.mark.parametrize("dim,prefix,geo_path", cases)
def test_timing(init_feelpp,dim,prefix,geo_path):
    fppc.Environment.changeRepository(
        directory=f"pyfeelpp-tests/timing/{prefix}")
    tic()
    print(f"geo {dim}D, file: {geo_path}")
    assert os.path.exists(geo_path)
    m = fppc.load(fppc.mesh(dim=dim,realdim=dim), geo_path, 0.1)
    toc(f"geo {dim}D")
    fppc.Environment.saveTimers(display=True)