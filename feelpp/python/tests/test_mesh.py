import feelpp.core as fppc
import sys
import pytest

@pytest.mark.parametrize("dim,order,realdim", [
    (1, 1, 1), (1, 2, 1), (1, 1, 2),
    (2, 1, 2), (2, 2, 2), (2, 1, 3),
    (3, 1, 3), (3, 2, 3),
])
def test_mesh_create(init_feelpp, dim, order, realdim):
    mesh_type = getattr(fppc, f"Mesh_S{dim}DG{order}R{realdim}")
    mesh = mesh_type.create(fppc.Environment.worldCommPtr())
    assert isinstance(mesh, mesh_type)
    assert mesh.dimension() == dim
    assert mesh.order() == order
    assert mesh.realDimension() == realdim


def run( m, geofile):
    if fppc.Environment.isMasterRank():
        print("mesh dim:", m.dimension())

    m = fppc.load(m, geofile, 0.1)
    if fppc.Environment.isMasterRank():
        print("mesh ", m.dimension(), "D nelts:", m.numGlobalElements())
        print("mesh ", m.dimension(), "D nfaces:", m.numGlobalFaces())
        print("mesh ", m.dimension(), "D hmin:", m.hMin())
        print("mesh ", m.dimension(), "D havg:", m.hAverage())
        print("mesh ", m.dimension(), "D hmax:", m.hMax())
        print("mesh ", m.dimension(), "D measure:", m.measure())


    
    r = fppc.elements(m)
    ne = fppc.nelements(r, True)
    if fppc.Environment.isMasterRank():
        print("mesh elts:", ne)
    r = fppc.boundaryfaces(m)
    nf = fppc.nelements(r, True)
    if fppc.Environment.isMasterRank():
        print("mesh boundary faces:", nf )
    s=fppc.createSubmesh(m,fppc.elements(m))  
    nes = fppc.nelements(fppc.elements(s), True)
    assert(ne==nes)
    s2 = fppc.createSubmesh(m, fppc.boundaryfaces(m))
    nfs = fppc.nelements(fppc.elements(s2), True)
    assert(nf == nfs)

cases = [
         (2,'feelpp2d','cases/feelpp2d/feelpp2d.geo'),
         (3,'feelpp3d','cases/feelpp2d/feelpp3d.geo')
        ]


#@pytest.mark.mpi
@pytest.mark.parametrize("dim,prefix,geo_path", cases)
def test_mesh(init_feelpp,dim,prefix,geo_path):
    fppc.Environment.changeRepository(directory=f"pyfeelpp-tests/mesh/{prefix}")
    print(geo_path)
    run( fppc.mesh(dim=dim,realdim=dim), geo_path )


