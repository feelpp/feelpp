import feelpp
import sys
e=feelpp.Environment(sys.argv)

geo={
    '2':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0],
    '3':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp3d/feelpp3d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0]
}

def run( m, geo ):
    m2d = feelpp.load(m,geo,0.1)
    Xh=feelpp.functionSpace(mesh=m2d)

    if e.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3=Xh.mesh()

    assert m3==m2d

    u=Xh.element()

    u.on(range=feelpp.elements(m2d),expr=feelpp.expr("x:x"))
    w=(2*u-u)-u
    l2_w=feelpp.normL2(range=feelpp.elements(m2d),expr=w)
    print("norm(w)",l2_w)
    assert abs(l2_w) < 1e-12
    assert u.functionSpace() == Xh
    assert u.size() == Xh.nDof()

run( feelpp.mesh(dim=2), geo['2'] )
run( feelpp.mesh(dim=3,realdim=3), geo['3'] )
