import feelpp
import sys
import pytest
from feelpp.timing import tic,toc

#@pytest.mark.mpi
def test_timing(init_feelpp):
    geo={
        '2':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0],
        '3':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp3d/feelpp3d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0]
    }
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/timing/")
    tic()
    g=geo['2']
    toc("geo 2D")
    tic()
    g=geo['3']
    toc("geo 3D")