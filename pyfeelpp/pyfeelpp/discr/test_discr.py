from pyfeelpp import core
import sys

e=core.Environment(sys.argv)

from pyfeelpp import mesh,discr,vf
geo={
    '2':core.download( "github:{repo:feelpp,path:quickstart/laplacian/feelpp2d/feelpp2d.geo}", worldComm=core.Environment.worldCommPtr() )[0],
    '3':core.download( "github:{repo:feelpp,path:quickstart/laplacian/feelpp3d/feelpp3d.geo}", worldComm=core.Environment.worldCommPtr() )[0]
}

def run( m, geo ):
    m2d = mesh.load(m,geo,0.1)
    Xh=discr.functionSpace(mesh=m2d)

    if e.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3=Xh.mesh()

    assert m3==m2d

    u=Xh.element()
    u.on(range=mesh.elements(m2d),expr=vf.expr("x:x"))

    assert u.functionSpace() == Xh
    assert u.size() == Xh.nDof()

run( mesh.mesh(dim=2), geo['2'] )
run( mesh.mesh(dim=3,realdim=3), geo['3'] )
