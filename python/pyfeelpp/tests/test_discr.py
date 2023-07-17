import feelpp
import sys
import pytest
from feelpp.timing import tic,toc


def run(m, geo):
    tic()
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    m2d= feelpp.load(m, mesh_name, 0.1)
    toc("Mesh")
    tic()
    Xh = feelpp.functionSpace(mesh=m2d)
    toc("functionSpace")

    if feelpp.Environment.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    assert m3 == m2d

    tic()
    u = Xh.element()
    toc("FunctionSpace::Element")
    tic()
    u.on(range=feelpp.elements(m2d), expr=feelpp.expr("x:x"))
    toc("FunctionSpace::Element::on")
    tic()
    w = (2*u-u)-u
    l2_w = feelpp.normL2(range=feelpp.elements(m2d), expr=w)
    toc("FunctionSpace::Element::normL2")
    print("norm(w)", l2_w)
    assert abs(l2_w) < 1e-12
    assert u.functionSpace() == Xh
    assert u.size() == Xh.nDof()

def run_vectorial(m, geo):
    tic()
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    m2d= feelpp.load(m, mesh_name, 0.1)
    toc("Mesh")
    tic()
    Xh = feelpp.functionSpace(mesh=m2d,space="Pchv")
    toc("functionSpace")

    if feelpp.Environment.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    assert m3 == m2d

    tic()
    u = Xh.element()
    toc("FunctionSpace::Element")
    tic()
    u.on(range=feelpp.elements(m2d), expr=feelpp.expr("{x,y}:x:y" if dim==2 else "{x,y,z}:x:y:z",row=dim,col=1))
    v=feelpp.vonmises(displacement=u,model={ "model": { "type": "linear-elasticity", "parameters": { "mu": 5.55555, "lambda": 12.96 } } } )
    toc("FunctionSpace::Element::on")
    tic()
    w = (2*u-u)-u
    l2_w = feelpp.normL2(range=feelpp.elements(m2d), expr=w)
    toc("FunctionSpace::Element::normL2")
    print("norm(w)", l2_w)
    assert abs(l2_w) < 1e-12
    assert u.functionSpace() == Xh
    assert u.size() == Xh.nDof()

def run_element(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    
    m2d= feelpp.load(m, mesh_name, 0.1)

    Xh = feelpp.functionSpace(mesh=m2d)

    if feelpp.Environment.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    assert m3 == m2d

    b = feelpp.backend(worldcomm=feelpp.Environment.worldCommPtr())
    v = b.newVector(dm=Xh.mapPtr())
    
    # test constant
    v.setConstant(1.0)
    # create u which use the same memory storage as v
    u = Xh.element(v)
    w = (2*u-u)-u
    l2_w = feelpp.normL2(range=feelpp.elements(m2d), expr=w)
    print("norm(w)", l2_w)
    assert abs(l2_w) < 1e-12
    
    uu = Xh.element()
    uu.on(range=feelpp.elements(m2d), expr=feelpp.expr("x:x"))
    vv = uu.to_petsc()
    u.on(range=feelpp.elements(m2d), expr=feelpp.expr("x:x"))# now v contains also x
    w=uu-u
    l2_w = feelpp.normL2(range=feelpp.elements(m2d), expr=w)
    print("norm(uu-u)={}, uu.l2={}, u.l2={}".format(l2_w, uu.l2Norm(), u.l2Norm()))
    assert abs(l2_w) < 1e-12
    
    vv-=v
    l2_w = w.l2Norm()
    print("norm(vv_petsc-v_petsc)={}, vv.l2={}, v.l2={}".format(l2_w, vv.l2Norm(), v.l2Norm()))
    assert abs(l2_w) < 1e-12

def run_element_alg(m, geo):
    mesh_name, dim, e_meas, e_s_1, e_s_2, e_s_bdy=geo
    
    m2d = feelpp.load(m, mesh_name, 0.1)

    Xh = feelpp.functionSpace(mesh=m2d)

    if feelpp.Environment.isMasterRank():
        print("Xh basisname: ", Xh.basisName())
        print("Xh nDof: ", Xh.nDof())
        print("Xh nLocalDof: ", Xh.nLocalDof())
        print("Xh nLocalDofWithGhost: ", Xh.nLocalDofWithGhost())
        print("Xh nLocalDofWithoutGhost: ", Xh.nLocalDofWithoutGhost())

    m3 = Xh.mesh()

    assert m3 == m2d

    b = feelpp.backend(worldcomm=feelpp.Environment.worldCommPtr())
    v = b.newVector(dm=Xh.mapPtr())
    
    # test constant
    v.setConstant(1.0)
    # create u which use the same memory storage as v
    u = Xh.element(v)

    us = Xh.newVectors(3)
    for ui in us:
        ui.setConstant(1.0)
    
    r = u.mDot(us)
    assert r.shape == (3,)
    for i in range(3):
        assert( r[i] == us[i].sum() )


geo_cases=[(2, feelpp.create_rectangle),
            (3, feelpp.create_box)]
 
@pytest.mark.parametrize("dim,geo", geo_cases)
def test_discr(dim,geo,init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_{}d".format(dim))
    tic()
    run( feelpp.mesh(dim=dim,realdim=dim), geo() )
    toc(f"run {dim}D")

@pytest.mark.parametrize("dim,geo", geo_cases)
def test_discrv(dim,geo,init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_{}d".format(dim))
    tic()
    run_vectorial( feelpp.mesh(dim=dim,realdim=dim), geo() )
    toc(f"run {dim}D")

@pytest.mark.parametrize("dim,geo", geo_cases)
def test_element(dim,geo,init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_{}d_element".format(dim))
    run_element( feelpp.mesh(dim=dim, realdim=dim), geo(filename="boxelement" if dim==3 else "rectelement") )

@pytest.mark.parametrize("dim,geo", geo_cases)
def test_element_alg(dim, geo, init_feelpp):
    feelpp.Environment.changeRepository(
        directory="pyfeelpp-tests/discr/test_{}d_element".format(dim))
    run_element_alg( feelpp.mesh(dim=dim, realdim=dim), geo(filename="boxelement" if dim==3 else "rectelement") )
