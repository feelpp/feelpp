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
    nf = feelpp.nelements(r, True)
    if feelpp.Environment.isMasterRank():
        print("mesh boundary faces:", nf )
    s=feelpp.createSubmesh(m,feelpp.elements(m))  
    nes = feelpp.nelements(feelpp.elements(s), True)
    assert(ne==nes)
    s2 = feelpp.createSubmesh(m, feelpp.boundaryfaces(m))
    nfs = feelpp.nelements(feelpp.elements(s2), True)
    assert(nf == nfs)

cases = [
         (2,'feelpp2d','cases/feelpp2d/feelpp2d.geo'),
         (3,'feelpp3d','cases/feelpp2d/feelpp3d.geo')
        ]


#@pytest.mark.mpi
@pytest.mark.parametrize("dim,prefix,geo_path", cases)
def test_mesh(init_feelpp,dim,prefix,geo_path):
    feelpp.Environment.changeRepository(directory=f"pyfeelpp-tests/mesh/{prefix}")
    print(geo_path)
    run( feelpp.mesh(dim=dim,realdim=dim), geo_path )


