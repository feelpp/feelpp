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
    ne = feelpp.nelements(r, True)
    if feelpp.Environment.isMasterRank():
        print("mesh elts:", ne)
    r = feelpp.boundaryfaces(m)
    nf = feelpp.nfaces(r, True)
    if feelpp.Environment.isMasterRank():
        print("mesh boundary faces:", nf )
    s=feelpp.createSubmesh(m,feelpp.elements(m))  
    nes = feelpp.nelements(feelpp.elements(s), True)
    assert(ne==nes)
#    s2 = feelpp.createSubmesh(m, feelpp.boundaryfaces(m))
#    nfs = feelpp.nfaces(feelpp.elements(s2), True)
#    assert(nf == nfs)

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

