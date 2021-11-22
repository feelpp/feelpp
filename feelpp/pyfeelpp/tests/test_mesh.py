import feelpp
import sys
import pytest

def run( m, geofile):
    if feelpp.Environment.isMasterRank():
        print("mesh dim:", m.dimension())

    m = feelpp.load(m, geofile, 0.1)
    if feelpp.Environment.isMasterRank():
        print("mesh ", m.dimension(), "D nelts:", m.numGlobalElements())
        print("mesh ", m.dimension(), "D nfaces:", m.numGlobalFaces())
        print("mesh ", m.dimension(), "D hmin:", m.hMin())
        print("mesh ", m.dimension(), "D havg:", m.hAverage())
        print("mesh ", m.dimension(), "D hmax:", m.hMax())
        print("mesh ", m.dimension(), "D measure:", m.measure())

    r = feelpp.elements(m)
    print("mesh elts:", feelpp.nelements(r, True))
    r = feelpp.boundaryfaces(m)
    print("mesh boundary faces:", feelpp.nfaces(r, True))
        

#@pytest.mark.mpi
def test_mesh(init_feelpp):
    geo={
        '2':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0],
        '3':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp3d/feelpp3d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0]
    }
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/mesh/test_2d")
    run( feelpp.mesh(dim=2), geo['2'] )
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/mesh/test_3d")
    run( feelpp.mesh(dim=3,realdim=3), geo['3'] )

