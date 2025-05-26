import re
import pytest
from feelpp.sympy.systems_pde import generate_stokes, generate_darcy

def test_stokes_2d():
    out = generate_stokes(nu=2, u_expr="{x,y}", p_expr="x+y", dim=2)
    # f should be a 2-vector literal with rank1
    f = out['f'].replace(' ', '')
    assert f.startswith('{')
    assert f.endswith(':rank1')
    # g = div(u) = div({x,y}) = 1+1 = 2
    assert out['g'].startswith('2'), f"Expected g starting with '2', got {out['g']}"
    # nu preserved
    assert out['nu'].startswith('2'), f"Expected nu='2', got {out['nu']}"

def test_darcy_2d():
    out = generate_darcy(kappa=3, p_expr="x**2+y", dim=2)
    # u = -3*[2*x,1] = {-6*x,-3}
    u = out['u'].replace(' ', '')
    assert u.startswith('{-6*x,-3}'), f"Unexpected u: {out['u']}"
    # q = div(u) = ∂x(-6*x) + ∂y(-3) = -6 + 0 = -6
    assert out['q'].startswith('-6'), f"Unexpected q: {out['q']}"
    # kappa preserved
    assert out['kappa'].startswith('3'), f"Expected kappa='3', got {out['kappa']}"