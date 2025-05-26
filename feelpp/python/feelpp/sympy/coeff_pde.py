#!/usr/bin/env python3
# feelpp/sympy/coeff_pde.py

import argparse
from sympy import symbols, sympify, Array
from .utils     import parse_vector, dot
from .operators import grad, div, syms, nsyms, n, dt, d2t
from .codegen   import sympytoginac


def generate_coeff_pde(
    d       = 0,
    d2      = 0,
    c       = 1,
    alpha   = "{0,0}",
    gamma   = "{0,0}",
    beta    = "{0,0}",
    a       = 0,
    p_expr  = None,
    dim     = 2,
    include_time = True
):
    """
    Build GiNaC/C99 code for the coefficient-form PDE,
    supporting both scalar and vector u.

    Generic form:
      d2 ∂²ₜu + d ∂ₜu
      + ∇·(−c ∇u - α u + γ)
      + β·∇u + a u = f
    """
    # 1) Base symbols
    x, y, z, t = symbols('x y z t')

    # 2) Parse u (scalar expr or vector literal)
    if isinstance(p_expr, str) and p_expr.strip().startswith('{'):
        u = parse_vector(p_expr, {'x':x,'y':y,'z':z,'t':t})
    else:
        u = sympify(p_expr, {'x':x,'y':y,'z':z,'t':t}) if p_expr else x + y

    # context for sympifying coefficients
    locals_map = {'x':x,'y':y,'z':z,'t':t,'u':u}

    # 3) Sympify scalars
    d  = sympify(d,  locals_map)
    d2 = sympify(d2, locals_map)
    c  = sympify(c,  locals_map)
    a  = sympify(a,  locals_map)

    # 4) Parse vectors (keep originals for codegen)
    alpha0 = parse_vector(alpha, locals_map)
    gamma0 = parse_vector(gamma, locals_map)
    beta0  = parse_vector(beta,  locals_map)

    # use working copies
    alpha = alpha0
    gamma = gamma0
    beta  = beta0

    # 5) Dimension checks
    for vec,name in ((alpha,'alpha'), (gamma,'gamma'), (beta,'beta')):
        if len(vec) != dim:
            raise ValueError(f"{name}-vector length {len(vec)} ≠ dim={dim}")

    # 6) Spatial operators
    S  = syms(dim)
    NS = nsyms(dim)
    grad_u    = grad(u, S)
    diff_flux = -c * grad_u

    # 7) Conservative convective flux: elementwise for vector
    if isinstance(u, Array):
        cons_flux = Array([alpha[i]*u[i] for i in range(dim)]) * -1
    else:
        cons_flux = -alpha * u

    # 8) Promote cons_flux & gamma to rank-2 if diff_flux is rank-2
    if isinstance(diff_flux, Array) and diff_flux.rank() == 2:
        if isinstance(cons_flux, Array) and cons_flux.rank() == 1:
            cons_flux = Array([cons_flux for _ in range(dim)])
        if isinstance(gamma, Array) and gamma.rank() == 1:
            gamma = Array([gamma for _ in range(dim)])

    # 9) Total conservative flux
    Fc = diff_flux + cons_flux + gamma

    # 10) Advective term
    if isinstance(u, Array):
        adv_term = Array([dot(beta, grad_u[i]) for i in range(dim)])
    else:
        adv_term = dot(beta, grad_u)

    # 11) Time terms
    if isinstance(u, Array):
        # vector case: compute each component
        if include_time:
            dt_u  = dt(u)
            d2t_u = d2t(u)
            term1 = Array([d*dt_u[i]   for i in range(dim)])
            term2 = Array([d2*d2t_u[i] for i in range(dim)])
        else:
            term1 = Array([0]*dim)
            term2 = Array([0]*dim)
    else:
        # scalar case
        if include_time:
            term1 = d * dt(u)
            term2 = d2 * d2t(u)
        else:
            term1 = sympify(0)
            term2 = sympify(0)

    # 12) Build RHS and normal flux
    divFc = div(Fc, S)
    if isinstance(u, Array):
        f_rhs = Array([
            term2[i] + term1[i] + divFc[i] + adv_term[i] + a*u[i]
            for i in range(dim)
        ])
    else:
        f_rhs = term2 + term1 + divFc + adv_term + a*u
    un = n(Fc, 1, NS)

    # 13) Emit GiNaC/C99 literals (using original vects)
    return {
        'd':     sympytoginac(d),
        'd2':    sympytoginac(d2),
        'c':     sympytoginac(c),
        'alpha': sympytoginac(alpha0),
        'gamma': sympytoginac(gamma0),
        'beta':  sympytoginac(beta0),
        'a':     sympytoginac(a),
        'f':     sympytoginac(f_rhs),
        'un':    sympytoginac(un)
    }

# -----------------------------------------------------------------------------
# Convenience wrappers
# -----------------------------------------------------------------------------

def generate_laplacian(c=1, p_expr=None, dim=2):
    return generate_coeff_pde(c=c, p_expr=p_expr, dim=dim, include_time=False)

def generate_convection_diffusion(c=1, beta="{1,0}", p_expr=None, dim=2):
    return generate_coeff_pde(c=c, beta=beta, p_expr=p_expr, dim=dim, include_time=False)

def generate_adr(c=1, beta="{1,0}", alpha="{0,0}", a=0, p_expr=None, dim=2):
    return generate_coeff_pde(c=c, beta=beta, alpha=alpha, a=a,
                              p_expr=p_expr, dim=dim, include_time=False)

def generate_wave(d=0, d2=1, c=1, alpha="{0,0}", gamma="{0,0}",
                  beta="{0,0}", a=0, p_expr=None, dim=2):
    return generate_coeff_pde(d=d, d2=d2, c=c,
                              alpha=alpha, gamma=gamma,
                              beta=beta, a=a,
                              p_expr=p_expr, dim=dim,
                              include_time=True)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate coefficient-form PDE code for Feel++"
    )
    subs = parser.add_subparsers(dest="cmd", required=True)

    # generic
    gen = subs.add_parser("coeff", help="generic PDE")
    gen.add_argument("--d",   type=str, default="0",   help="∂ₜu coeff")
    gen.add_argument("--d2",  type=str, default="0",   help="∂²ₜu coeff")
    gen.add_argument("--c",   type=str, default="1",   help="diffusion c")
    gen.add_argument("--alpha", type=str, default="{0,0}", help="vector α")
    gen.add_argument("--gamma", type=str, default="{0,0}", help="vector γ")
    gen.add_argument("--beta",  type=str, default="{0,0}", help="vector β")
    gen.add_argument("--a",     type=str, default="0",   help="reaction a")
    gen.add_argument("--p",     type=str, default=None,  help="u(x,t) expr")
    gen.add_argument("--dim",   type=int, default=2,     help="space dim")

    # laplacian
    lap = subs.add_parser("laplacian", help="pure Laplace")
    lap.add_argument("--c",   type=str, default="1",    help="diffusion c")
    lap.add_argument("--p",   type=str, default=None,   help="u(x,t) expr")
    lap.add_argument("--dim", type=int, default=2,      help="space dim")

    # convection-diffusion
    cd = subs.add_parser("cd", help="convection–diffusion")
    cd.add_argument("--c",   type=str, default="1",      help="diffusion c")
    cd.add_argument("--beta",type=str, default="{1,0}",  help="vector β")
    cd.add_argument("--p",   type=str, default=None,     help="u(x,t) expr")
    cd.add_argument("--dim", type=int, default=2,        help="space dim")

    # ADR
    adr = subs.add_parser("adr", help="advection-diffusion-reaction")
    adr.add_argument("--c",   type=str, default="1",     help="diffusion c")
    adr.add_argument("--beta",type=str, default="{1,0}", help="vector β")
    adr.add_argument("--alpha",type=str,default="{0,0}", help="vector α")
    adr.add_argument("--a",   type=str, default="0",     help="reaction a")
    adr.add_argument("--p",   type=str, default=None,    help="u(x,t) expr")
    adr.add_argument("--dim", type=int, default=2,       help="space dim")

    # wave
    wave = subs.add_parser("wave", help="damped wave")
    wave.add_argument("--d",   type=str, default="0",    help="∂ₜu coeff")
    wave.add_argument("--d2",  type=str, default="1",    help="∂²ₜu coeff")
    wave.add_argument("--c",   type=str, default="1",    help="diffusion c")
    wave.add_argument("--alpha",type=str,default="{0,0}",help="vector α")
    wave.add_argument("--gamma",type=str,default="{0,0}",help="vector γ")
    wave.add_argument("--beta", type=str, default="{0,0}",help="vector β")
    wave.add_argument("--a",   type=str, default="0",    help="reaction a")
    wave.add_argument("--p",   type=str, default=None,   help="u(x,t) expr")
    wave.add_argument("--dim", type=int, default=2,      help="space dim")

    args = parser.parse_args(argv)

    if args.cmd == "coeff":
        out = generate_coeff_pde(
            d=args.d, d2=args.d2,
            c=args.c, alpha=args.alpha, gamma=args.gamma,
            beta=args.beta, a=args.a,
            p_expr=args.p, dim=args.dim,
            include_time=True
        )
    elif args.cmd == "laplacian":
        out = generate_laplacian(c=args.c, p_expr=args.p, dim=args.dim)
    elif args.cmd == "cd":
        out = generate_convection_diffusion(
            c=args.c, beta=args.beta,
            p_expr=args.p, dim=args.dim
        )
    elif args.cmd == "adr":
        out = generate_adr(
            c=args.c, beta=args.beta,
            alpha=args.alpha, a=args.a,
            p_expr=args.p, dim=args.dim
        )
    else:  # wave
        out = generate_wave(
            d=args.d, d2=args.d2,
            c=args.c, alpha=args.alpha, gamma=args.gamma,
            beta=args.beta, a=args.a,
            p_expr=args.p, dim=args.dim
        )

    for name, val in out.items():
        print(f"{name}: {val}")

if __name__ == "__main__":
    main()