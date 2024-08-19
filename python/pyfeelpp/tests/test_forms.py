import feelpp.core as fppc
import sys
import pytest
from feelpp.core.timing import tic,toc
from feelpp.core.forms import *
import numpy as np

def run(m, geo):
    tic()
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    m2d= fppc.load(m, mesh_name, 0.1)
    toc("Mesh")
    tic()
    Xh = fppc.functionSpace(mesh=m2d)
    toc("functionSpace")

    if fppc.Environment.isMasterRank():
        print("[test_forms] Xh basisname: ", Xh.basisName())
        print("[test_forms] Xh nDof: ", Xh.nDof())
        print("[test_forms] Xh nLocalDof: ", Xh.nLocalDof())
        print("[test_forms] Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("[test_forms] Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    a = mass(test=Xh,trial=Xh)
    u=Xh.element()
    u.setConstant(1)

    al = a(u)
    print(f"[form1] measure using l={al(u)}")
    print(f"[form2] measure={a(u,u)}")
    assert np.isclose(a(u,u),e_meas)

    l = flux(test=Xh)
    f=Xh.element()
    f.setConstant(1)
    print(f"[form1] measure using l={l(f)}")
    assert np.isclose(l(f),e_meas)
    l_bdy=flux(test=Xh,boundary=True)
    print(f"[form1] measure using l_bdy={l_bdy(f)}")
    assert np.isclose(l_bdy(f),e_s_bdy)

    b=aGradGrad(test=Xh,trial=Xh)
    l.zero()
    b.on(rhs=l,expr="x:x",element=u)
    b.solve(rhs=l,solution=u)
    # u should be equal to 1
    v=Xh.element()
    v.on(range=fppc.elements(Xh.mesh()), expr=fppc.expr("x:x:y"))
    l2_e = fppc.normL2(range=fppc.elements(Xh.mesh()), expr=u-v)
    print(f"norm(u-v)={l2_e}") # should be 0 up to machine precision
    assert np.isclose(l2_e,0)

geo_cases=[(2, fppc.create_rectangle),
            (3, fppc.create_box)]
 
@pytest.mark.parametrize("dim,geo", geo_cases)
def test_form(dim,geo,init_feelpp):
    fppc.Environment.changeRepository(
        directory="pyfeelpp-tests/forms/test_{}d".format(dim))
    tic()
    run( fppc.mesh(dim=dim,realdim=dim), geo() )
    toc(f"run {dim}D")


