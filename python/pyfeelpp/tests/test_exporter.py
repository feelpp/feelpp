import sys
import feelpp.core as fppc 
import pytest

def run(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy = geo()
    m2d = fppc.load(m, mesh_name, 0.1)

    Xh = fppc.functionSpace(space="Pch", mesh=m2d, order=1)
    Xhv = fppc.functionSpace(space="Pchv", mesh=m2d, order=1)
    P0h = fppc.functionSpace(space="Pdh", mesh=m2d, order=0)
    #u=Xh.elementFromExpr("{sin(2*pi*x)*cos(pi*y)}:x:y")
    u = Xh.element()
    u.on(range=fppc.elements(m2d), expr=fppc.expr("x*x:x"))
    v = Xhv.element()
    v.on(range=fppc.elements(m2d), expr=fppc.expr("{x,y}:x:y" if dim==2 else "{x,y,z}:x:y:z",row=dim))

    e = fppc.exporter(mesh=m2d, name="feelpp"+str(m.dimension())+"d")
    e.add("un", 1.)
    e.add("u", u)
    e.add("v", v)
    e.add("pid", fppc.pid(P0h))
    e.save()


def test_exporter(init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/exporter")
    geo={
        '2':fppc.create_rectangle,
        '3':fppc.create_box
    }
    run( fppc.mesh( dim=2 ), geo['2'] )
    run( fppc.mesh( dim=3, realdim=3 ), geo['3'] )



