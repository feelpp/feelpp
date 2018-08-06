import sys
from pyfeelpp import *
e=core.Environment(sys.argv)
import pyfeelpp.discr as discr
import pyfeelpp.filters as filters


m2d=mesh.Mesh_2()
m2d = mesh.load(m2d,"triangle.geo",0.1)

Xh=discr.Pch_2D_P1(mesh=m2d)
P0h = discr.Pdh_2D_P0(mesh=m2d)
u=Xh.elementFromExpr("{sin(2*pi*x)*cos(pi*y)}:x:y")
e = filters.exporter(mesh=m2d)
e.addScalar("un", 1.)
e.addP1c("u",u);
e.addP0d("pid",discr.pid( P0h ));
e.save()
