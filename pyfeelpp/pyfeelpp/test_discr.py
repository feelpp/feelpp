import core
import mesh
import sys,time
e=core.Environment(sys.argv)

import discr



print("pid:",e.worldComm().localRank() )
m2d=mesh.Mesh_2()
m2d = mesh.load(m2d,"triangle.geo",0.1)

Xh=discr.Pch_2D_P1(mesh=m2d)
Yh=discr.Pch_2D_P1(mesh=m2d)


print("Xh basisname: ", Xh.basisName())
print("Xh nDof: ", Xh.nDof())
print("Xh nLocalDof: ", Xh.nLocalDof())
print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

m3=Xh.mesh()

assert m3==m2d

u=Xh.element()

assert u.functionSpace() == Xh
assert u.size() == Xh.nDof()
assert u.functionSpace() != Yh

#help(u)
