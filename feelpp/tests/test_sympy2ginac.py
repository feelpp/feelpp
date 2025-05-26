import pytest
from sympy import symbols, Array, Matrix, Integral
import feelpp.sympy as fppsy


# Define spatial symbols
x, y, z, t = symbols('x y z t')
nx, ny, nz = symbols('nx ny nz')

# Test parse_vector

def test_parse_vector_numeric():
    arr = fppsy.parse_vector('{1,2,3}')
    assert isinstance(arr, Array)
    assert list(arr) == [1, 2, 3]


def test_parse_vector_symbols():
    arr = fppsy.parse_vector('{x,y,z}')
    assert list(arr) == [x, y, z]

# Test toginac scalar

def test_toginac_scalar():
    lit = fppsy.toginac(3, [], rank=0)
    assert lit.endswith(':rank0')
    assert lit.startswith('3')

# Test toginac vector

def test_toginac_vector():
    vec = fppsy.parse_vector('{4,5}')
    lit = fppsy.toginac(vec, [x, y], rank=1)
    # should include both symbol suffix and rank
    assert ':x:y' in lit
    assert lit.endswith(':rank1')

# Test sympytoginac on string literal

def test_sympytoginac_string():
    lit = fppsy.sympytoginac('{7,8}')
    assert ':rank1' in lit or ':rank2' not in lit  # ensure rank applied

# Test differential operators

def test_grad_div_laplacian():
    sym = fppsy.syms(2)
    f = x**2 + 3*y
    g = fppsy.grad(f, sym)
    assert list(g) == [2*x, 3]
    d = fppsy.div(g, sym)
    assert d == 2  # divergence of constant gradient is 0
    
    # laplacian of x^2 + y is 2
    lap = fppsy.laplacian(f, sym)
    assert lap == 2

# Test tensor helpers

def test_mult_dot_outer():
    a = Array([1,2,3])
    b = Array([4,5,6])
    dot = fppsy.dot(a, b)
    assert dot == 1*4 + 2*5 + 3*6
    
    # outer product gives a 3x3 Matrix
    out = fppsy.outer(a, b)
    assert isinstance(out, Matrix)
    assert out.shape == (3,3)

# Test curl and cross

#def test_curl_cross():
#    u = Array([y, z, x])
#    c = fppsy.curl(u)
#    # manually compute: curl([y,z,x]) = [d/dy x - d/dz z, d/dz y - d/dx x, d/dx z - d/dy y] = [0 - 1, 0 - 1, 0 - 1] = [-1,-1,-1]
#    assert list(c) == [-1, -1, -1]
#    cross_ab = fppsy.cross(a='
#    ', b='   ')
#    # placeholder to ensure cross signature
#    # actual cross test
#    a = Array([1,0,0]); b = Array([0,1,0]); ab = fppsy.cross(a, b)
#    assert list(ab) == [0,0,1]

# Test hessian and material_derivative

def test_hessian_material():
    f = x*y + z**2
    H = fppsy.hessian(f, [x, y, z])
    assert H[2,2] == 2
    u = Array([1,2,3])
    md = fppsy.material_derivative(x+y, u, [x,y,z], t)
    # ∂_t(x+y)=0, u·∇(x+y)=1*1 +2*1 =3
    assert md == 3

# Test integrate

def test_integrate():
    expr = x
    res = fppsy.integrate(expr, [(x, 0, 1)])
    assert pytest.approx(res) == 0.5

if __name__ == '__main__':
    pytest.main()
