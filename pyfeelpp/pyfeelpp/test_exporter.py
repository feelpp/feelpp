import core
import mesh
import sys,time
e=core.Environment(sys.argv)

import discr
import exporter


m2d=mesh.Mesh_2()
m2d = mesh.load(m2d,"triangle.geo",0.1)

Xh=discr.Pch_2D_P1(mesh=m2d)
u=Xh.element()

e = exporter.exporter(mesh=m2d)
e.add()
e.save()
