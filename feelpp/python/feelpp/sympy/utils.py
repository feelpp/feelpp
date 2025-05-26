# feelpp/sympy/utils.py

import re
from sympy import sympify, Array, Integral

def get_var(name, default, cast=None, context=None):
    ctx = context if context is not None else locals()
    val = ctx.get(name, default)
    if val is None or val == '':
        val = default
    return cast(val) if (cast and val is not None) else val

def parse_vector(expr_str, locals_map=None):
    """
    Parse a Feel++ vector literal into a Sympy Array,
    sympifying each entry with optional locals_map.
    """
    s = str(expr_str).strip()
    m = re.match(r'^\{([^}]*)\}', s)
    if not m:
        raise ValueError(f"cannot parse vector from `{expr_str}`")
    body = m.group(1)
    parts = body.split(',') if ',' in body else body.split(':')
    if locals_map:
        return Array([sympify(p, locals_map) for p in parts])
    else:
        return Array([sympify(p) for p in parts])

def dot(a, b):
    return sum(a[i]*b[i] for i in range(len(a)))

def outer(a, b):
    from sympy import Matrix
    from sympy.tensor.array import tensorproduct
    return Matrix(tensorproduct(a, b))

def integrate(f, limits):
    res = f
    for sym, lo, hi in limits:
        res = Integral(res, (sym, lo, hi))
    return res.doit()