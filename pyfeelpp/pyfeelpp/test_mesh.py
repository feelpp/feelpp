import core
import mesh
import sys,time

e=core.Environment(sys.argv)
m1d=mesh.Mesh_1()
m2d=mesh.Mesh_2()
if e.isMasterRank():
    print("mesh dim:", m1d.dimension())
    print("mesh dim:", m2d.dimension())

m2d=mesh.load(m2d,"triangle.geo",2)
help(m2d)
if e.isMasterRank():
    print("mesh 2D nelts:", m2d.numGlobalElements() )
    print("mesh 2D nfaces:", m2d.numGlobalFaces() )
    print("mesh 2D hmin:", m2d.hMin())
    print("mesh 2D havg:", m2d.hAverage())
    print("mesh 2D hmax:", m2d.hMax())
    print("mesh 2D measure:", m2d.measure())
    
    
r = mesh.elements(m2d)
print("mesh elts:", mesh.nelements(r,True))
r = mesh.boundaryfaces(m2d)
print("mesh boundary faces:", mesh.nfaces(r,True))
