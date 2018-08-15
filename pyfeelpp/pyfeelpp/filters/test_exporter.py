import sys
from pyfeelpp import core

e=core.Environment(sys.argv)

from pyfeelpp import mesh,discr,filters,vf

rd = core.RemoteData("github:{repo:feelpp,path:quickstart/laplacian/feelpp2d/feelpp2d.geo}",worldComm=core.Environment.worldCommPtr())
if rd.canDownload():
    d=e.downloadsRepository()
    geo={'2':rd.download( d )[0]}
    if e.isMasterRank():
        print(geo)

rd = core.RemoteData("github:{repo:feelpp,path:quickstart/laplacian/feelpp3d/feelpp3d.geo}",worldComm=core.Environment.worldCommPtr())
if rd.canDownload():
    d=e.downloadsRepository()
    geo['3']=rd.download( d )[0]
    if e.isMasterRank():
        print(geo)        

def run( m, geo ):
    m2d = mesh.load(m,geo,0.1)

    Xh=discr.functionSpace( type="Pch", mesh=m2d, order=1 )
    P0h = discr.functionSpace( type="Pdh", mesh=m2d, order=0 )
    #u=Xh.elementFromExpr("{sin(2*pi*x)*cos(pi*y)}:x:y")
    u=Xh.element()
    u.on(range=mesh.elements(m2d),expr=vf.expr("x*x:x"))

    e = filters.exporter(mesh=m2d,name="feelpp"+str(m.dimension())+"d")
    e.addScalar("un", 1.)
    e.addP1c("u",u);
    e.addP0d("pid",discr.pid( P0h ));
    e.save()

run( mesh.mesh( dim=2 ), geo['2'] )
run( mesh.mesh( dim=3, realdim=3 ), geo['3'] )    
