# feelpp/sympy/operators.py
"""
Differential and tensor operators for Feel++ symbolic PDEs.
Supports both scalar Expr and Array fields, including first- and second-order time derivatives.
"""

from sympy import derive_by_array, transpose, diff, Matrix, Array
from sympy.tensor.array import tensorcontraction, tensorproduct
from .utils import dot
from . import x, y, z, t, nx, ny, nz

def dx(f):
    """∂f/∂x"""
    return derive_by_array(f, [x])

def dy(f):
    """∂f/∂y"""
    return derive_by_array(f, [y])

def dz(f):
    """∂f/∂z"""
    return derive_by_array(f, [z])

def dt(f):
    """
    First time derivative ∂ₜf.
    - If f is a scalar Expr, returns diff(f, t).
    - If f is an Array, returns an Array of componentwise diff.
    """
    from sympy import diff as _diff
    if isinstance(f, Array):
        return Array([_diff(f[i], t) for i in range(len(f))])
    else:
        return _diff(f, t)

def d2t(f):
    """
    Second time derivative ∂²ₜf.
    - If f is a scalar Expr, returns diff(f, t, 2).
    - If f is an Array, returns an Array of componentwise second derivatives.
    """
    from sympy import diff as _diff
    if isinstance(f, Array):
        return Array([_diff(f[i], t, 2) for i in range(len(f))])
    else:
        return _diff(f, t, 2)

def grad(f, syms_list):
    """Gradient ∇f over coordinates in syms_list"""
    return derive_by_array(f, syms_list)

def symgrad(f, syms_list):
    """Symmetric gradient ∇f + (∇f)ᵀ"""
    g = derive_by_array(f, syms_list)
    return g + transpose(g)

def div(f, syms_list):
    """Divergence of vector/tensor field f over syms_list"""
    return tensorcontraction(derive_by_array(f, syms_list), (0, 1))

def laplacian(f, syms_list):
    """Laplacian Δf = ∇·∇f"""
    return tensorcontraction(
        derive_by_array(derive_by_array(f, syms_list), syms_list),
        (0, 1)
    )

def curl(u):
    """Vector curl ∇×u for a 3D vector u"""
    return Array([
        dy(u[2]) - dz(u[1]),
        dz(u[0]) - dx(u[2]),
        dx(u[1]) - dy(u[0])
    ])

def hessian(f, syms_list=[x, y, z]):
    """Hessian matrix ∂²f/∂xi∂xj"""
    return Matrix([[diff(f, s1, s2) for s2 in syms_list] for s1 in syms_list])

def material_derivative(f, u, syms_list=[x, y, z], t_sym=t):
    """
    Material derivative D_t f = ∂ₜf + u·∇f.
    Works for scalar f and vector u.
    """
    term_t = diff(f, t_sym)
    grad_f = derive_by_array(f, syms_list)
    return term_t + dot(u, grad_f)

def syms(dim):
    """Return spatial symbols [x,y,z] for dimension dim"""
    return Array([x] if dim == 1 else [x, y] if dim == 2 else [x, y, z])

def nsyms(dim):
    """Return normal symbols [nx,ny,nz] for dimension dim"""
    return Array([nx] if dim == 1 else [nx, ny] if dim == 2 else [nx, ny, nz])

def n(a, c=1, ns=None):
    """
    Boundary normal flux n·(c a)
    For vector a, returns scalar; for Array a, returns Array contraction.
    """
    if ns is None:
        ns = nsyms(len(a))
    axis = 0 if Array(a).rank() == 1 else 1
    return tensorcontraction(tensorproduct(c * a, ns), (axis, axis + 1))