import sys
import feelpp 

e=feelpp.Environment(sys.argv)


geo={
    '2':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0],
    '3':feelpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp3d/feelpp3d.geo}", worldComm=feelpp.Environment.worldCommPtr() )[0]
}

def run( m, geo ):
    m2d = feelpp.load(m,geo,0.1)

    Xh=feelpp.functionSpace( space="Pch", mesh=m2d, order=1 )
    P0h = feelpp.functionSpace( space="Pdh", mesh=m2d, order=0 )
    #u=Xh.elementFromExpr("{sin(2*pi*x)*cos(pi*y)}:x:y")
    u=Xh.element()
    u.on(range=feelpp.elements(m2d),expr=feelpp.expr("x*x:x"))

    e = feelpp.exporter(mesh=m2d,name="feelpp"+str(m.dimension())+"d")
    e.addScalar("un", 1.)
    e.addP1c("u",u)
    e.addP0d("pid",feelpp.pid( P0h ))
    e.save()

run( feelpp.mesh( dim=2 ), geo['2'] )
run( feelpp.mesh( dim=3, realdim=3 ), geo['3'] )
