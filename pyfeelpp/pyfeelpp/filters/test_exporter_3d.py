import filters

e=core.Environment(sys.argv)

m=mesh.Mesh_3()
m = mesh.load(m,"feelpp3d.geo",0.1)

Xh=discr.Pch_3D_P1(mesh=m)
P0h = discr.Pdh_3D_P0(mesh=m)
u=Xh.elementFromExpr("{sin(2*pi*x)*cos(pi*y)*cos(pi*z)}:x:y:z")
e = exporter.exporter(mesh=m)
e.addScalar("un", 1.)
e.addP1c("u",u);
e.addP0d("pid",discr.pid( P0h ));
e.save()
