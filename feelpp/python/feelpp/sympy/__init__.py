# feelpp/sympy/__init__.py
"""
feelpp.sympy
=============
Symbolic PDE codegen utilities, wrapping Sympy → GiNaC/C99.
"""

from sympy import symbols

# Fundamental symbols for differential operators:
x, y, z, t = symbols("x y z t")
nx, ny, nz = symbols("nx ny nz")


from .utils      import get_var, parse_vector, dot, outer, integrate
from .codegen    import toginac, sympytoginac
from .operators  import grad, div, laplacian, curl, hessian, material_derivative, syms, nsyms, n
from .coeff_pde  import generate_coeff_pde, generate_laplacian, generate_convection_diffusion, generate_adr