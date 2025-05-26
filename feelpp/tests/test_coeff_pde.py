import re
import pytest
from feelpp.sympy.coeff_pde import (
    generate_laplacian,
    generate_convection_diffusion,
    generate_adr,
    generate_wave,
)

def assert_contains_all(s: str, substrs):
    """
    Assert that each substring is in s, ignoring whitespace differences.
    """
    # normalize by removing spaces
    s_norm = s.replace(" ", "")
    for sub in substrs:
        sub_norm = sub.replace(" ", "")
        assert sub_norm in s_norm, f"Expected '{sub}' in '{s}'"

def test_laplacian_scalar():
    out = generate_laplacian(c=1, p_expr="x+y", dim=2)
    # f should be zero
    assert out['f'].startswith("0")
    print(out['un'])
    # boundary flux un = -nx - ny
    assert_contains_all(out['un'], ["-nx", "-ny"])
    # ensure time coefficients are zero
    assert out['d'] == "0:rank0"
    assert out['d2'] == "0:rank0"

def test_convection_diffusion():
    out = generate_convection_diffusion(c=1, beta="{1,0}", p_expr="x+y", dim=2)
    # advective term = 1
    assert out['f'].startswith("1")
    # same un as laplacian
    assert_contains_all(out['un'], ["-nx", "-ny"])

def test_adr():
    out = generate_adr(c=1, beta="{1,0}", alpha="{0,0}", a=0, p_expr="x+y", dim=2)
    # identical to convection-diffusion in this simple case
    assert out['f'].startswith("1")
    assert_contains_all(out['un'], ["-nx", "-ny"])

def test_wave_scalar():
    out = generate_wave(
        d=0.5, d2=2.0, c=1,
        alpha="{0,0}", gamma="{0,0}",
        beta="{0,0}", a=0,
        p_expr="sin(t)*x", dim=2
    )
    # f should contain both sin(t) and cos(t) terms
    assert re.search(r"sin\(t\)", out['f'])
    assert re.search(r"cos\(t\)", out['f'])
    # un = -nx*sin(t)
    assert_contains_all(out['un'], ["-nx", "sin"])

def test_laplacian_vector_u():
    # vector solution u = {x, y}
    out = generate_laplacian(c=1, p_expr="{x,y}", dim=2)
    f = out['f']
    un = out['un']
    # f must be zero vector: {0,0}
    assert f.replace(' ', '').startswith('{0,0}'), f"Expected zero vector f, got {f}"
    # un must be vector of -nx, -ny and rank1
    un_norm = un.replace(' ', '')
    assert un_norm.startswith('{-nx,-ny}'), f"Unexpected un {un}"
    assert ':rank1' in un_norm, f"Expected rank1 for un, got {un}"

def test_adr_vector_u():
    # vector solution u = {x, y}, alpha zero, beta zero, a zero -> purely laplacian behavior
    out = generate_adr(c=1, beta="{0,0}", alpha="{0,0}", a="0", p_expr="{x,y}", dim=2)
    f = out['f']
    # should be zero vector
    assert f.replace(' ', '').startswith('{0,0}'), f"Expected zero vector f, got {f}"

def test_wave_vector_u_vector_solution():
    # wave with u vector, only second derivative term (d2=1)
    out = generate_wave(
        d="0", d2="1", c="1",
        alpha="{0,0}", gamma="{0,0}",
        beta="{0,0}", a="0",
        p_expr="{t,x}", dim=2
    )
    f = out['f']
    # f should be zero vector from second‐time derivative only
    assert f.replace(' ', '').startswith('{0,0}'), f"Expected zero vector f, got {f}"

    un = out['un'].replace(' ', '')
    # un should be a 2-vector with one zero entry and have rank1
    assert un.startswith('{') and ',0' in un, f"Unexpected un {un}"
    assert ':rank1' in un, f"Expected rank1 for un, got {un}"

if __name__ == "__main__":
    pytest.main()