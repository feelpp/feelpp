from pyfeelpp import core
import sys
e=core.Environment(sys.argv)
rd = core.RemoteData("github:{repo:feelpp,path:quickstart/laplacian/feelpp2d/feelpp2d.geo}",worldComm=core.Environment.worldCommPtr())
if rd.canDownload():
    d=e.downloadsRepository()
    print("download mesh in ", d);
    geo={'2':rd.download( d )[0]}
    print(geo)
rd = core.RemoteData("github:{repo:feelpp,path:quickstart/laplacian/feelpp3d/feelpp3d.geo}",worldComm=core.Environment.worldCommPtr())
if rd.canDownload():
    d=e.downloadsRepository()
    print("download mesh in ", d);
    geo['3']=rd.download( d )[0]
    print(geo)
from pyfeelpp import mesh

def run( m, geofile ):
    if e.isMasterRank():
        print("mesh dim:", m.dimension())
    
    m=mesh.load(m,geofile,0.1)
    if e.isMasterRank():
        print("mesh ",m.dimension(),"D nelts:", m.numGlobalElements() )
        print("mesh ",m.dimension(),"D nfaces:", m.numGlobalFaces() )
        print("mesh ",m.dimension(),"D hmin:", m.hMin())
        print("mesh ",m.dimension(),"D havg:", m.hAverage())
        print("mesh ",m.dimension(),"D hmax:", m.hMax())
        print("mesh ",m.dimension(),"D measure:", m.measure())
    
    r = mesh.elements(m)
    print("mesh elts:", mesh.nelements(r,True))
    r = mesh.boundaryfaces(m)
    print("mesh boundary faces:", mesh.nfaces(r,True))


run( mesh.mesh(dim=2), geo['2'] )
run( mesh.mesh(dim=3,realdim=3), geo['3'] )

