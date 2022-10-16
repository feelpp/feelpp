import feelpp
import sys
import pytest
from feelpp.timing import tic,toc


def run(m, geo):
    tic()
    m2d = feelpp.load(m, geo, 0.1)
    toc("Mesh")
    tic()
    Xh = feelpp.functionSpace(mesh=m2d)
    toc("functionSpace")

    if feelpp.Environment.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    assert m3 == m2d

    tic()
    u = Xh.element()
    toc("FunctionSpace::Element")
    tic()
    u.on(range=feelpp.elements(m2d), expr=feelpp.expr("x:x"))
    toc("FunctionSpace::Element::on")
    tic()
    w = (2*u-u)-u
    l2_w = feelpp.normL2(range=feelpp.elements(m2d), expr=w)
    toc("FunctionSpace::Element::normL2")
    print("norm(w)", l2_w)
    assert abs(l2_w) < 1e-12
    assert u.functionSpace() == Xh
    assert u.size() == Xh.nDof()


#@pytest.mark.mpi
def test_discr(init_feelpp):
    geo={
        '2':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0],
        '3':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp3d/feelpp3d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0]
    }
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_2d")
    tic()
    run( feelpp.mesh(dim=2), geo['2'] )
    toc("run 2D")
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_3d")
    tic()
    run( feelpp.mesh(dim=3,realdim=3), geo['3'] )
    toc("run 3D")
    
