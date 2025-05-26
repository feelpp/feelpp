#!/usr/bin/env python3
# feelpp/sympy/systems_pde.py

"""
Build GiNaC/C99 code for small PDE systems (Stokes, Darcy) via symbolic generation.
"""

from sympy import symbols, sympify, Array
from .utils     import parse_vector, dot
from .operators import grad, div, laplacian, syms
from .codegen   import sympytoginac

def generate_stokes(
    nu=1,               # viscosity
    u_expr=None,        # manufactured velocity vector "{u1,u2,...}"
    p_expr=None,        # manufactured pressure expression
    dim=2
):
    """
    Steady Stokes system:
      -nu Δu + ∇p = f   (vector)
       ∇·u       = g   (scalar)
    Returns dict with 'f', 'g', 'nu'.
    """
    # 1) symbols
    x, y, z = symbols('x y z')

    # 2) parse u and p
    default_u = '{' + ','.join('x' for _ in range(dim)) + '}'
    u = parse_vector(u_expr if u_expr is not None else default_u,
                     {'x':x, 'y':y, 'z':z})
    p = sympify(p_expr, {'x':x, 'y':y, 'z':z}) if p_expr is not None else x + y

    # 3) sympify nu with u,p in context
    nu_sym = sympify(nu, {'x':x,'y':y,'z':z,'u':u,'p':p})

    # 4) operators
    S = syms(dim)
    lap_u  = laplacian(u, S)    # vector Laplacian
    grad_p = grad(p, S)         # gradient of p

    # 5) momentum RHS and continuity RHS
    f_vec = -nu_sym * lap_u + grad_p
    g_sca = div(u, S)

    # 6) emit
    return {
        'f':   sympytoginac(f_vec),   # vector literal
        'g':   sympytoginac(g_sca),   # scalar literal
        'nu':  sympytoginac(nu_sym)
    }


def generate_darcy(
    kappa=1,           # permeability inverse
    p_expr=None,       # manufactured pressure
    dim=2
):
    """
    Darcy flow (steady):
      u = -kappa ∇p   (vector)
      div(u) = q      (scalar)
    Returns dict with 'u', 'q', 'kappa'.
    """
    # 1) symbols
    x, y, z = symbols('x y z')

    # 2) parse p
    p = sympify(p_expr, {'x':x,'y':y,'z':z}) if p_expr is not None else x + y

    # 3) sympify kappa
    kappa_sym = sympify(kappa, {'x':x,'y':y,'z':z,'p':p})

    # 4) operators
    S = syms(dim)
    grad_p = grad(p, S)
    u_vec  = -kappa_sym * grad_p
    q_sca  = div(u_vec, S)

    # 5) emit
    return {
        'u':     sympytoginac(u_vec),    # vector literal
        'q':     sympytoginac(q_sca),    # scalar literal
        'kappa': sympytoginac(kappa_sym)
    }