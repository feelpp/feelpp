# feelpp/sympy/api.py
"""
High-level Python API for generating PDE coefficients.
"""
from .coeff_pde import (
    generate_coeff_pde,
    generate_laplacian,
    generate_convection_diffusion,
    generate_adr,
    generate_wave,
)
from .systems_pde import generate_stokes, generate_darcy


def get_coefficients(pde: str, **kwargs) -> dict:
    """
    Fetch GiNaC/C99 literals for a given PDE form.

    Args:
      pde     : one of 'coeff', 'laplacian', 'cd', 'adr', 'wave', 'stokes', 'darcy'
      kwargs  : parameters for the generator (as in the underlying functions)

    Returns:
      dict mapping coefficient names to GiNaC literal strings.
    """
    if pde == "coeff":
        return generate_coeff_pde(**kwargs)
    elif pde == "laplacian":
        return generate_laplacian(**kwargs)
    elif pde == "cd":
        return generate_convection_diffusion(**kwargs)
    elif pde == "adr":
        return generate_adr(**kwargs)
    elif pde == "wave":
        return generate_wave(**kwargs)
    elif pde == "stokes":
        return generate_stokes(**kwargs)
    elif pde == "darcy":
        return generate_darcy(**kwargs)
    else:
        raise ValueError(f"Unknown PDE type '{pde}'")
