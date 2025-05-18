from sympy import *
import re

# -----------------------------------------------------------------------------
# Core utilities and code-generation for Feel++ PDE expressions via GiNaC/C99
# -----------------------------------------------------------------------------

# Fundamental symbols
x, y, z, t = symbols("x y z t")
nx, ny, nz = symbols("nx ny nz")
alpha, T = symbols("alpha T")

# ----------------------------------------------------------------------------
# Utility: fetch variable with default, and parse Feel++-style vector literals
# ----------------------------------------------------------------------------

def get_var(name, default, cast=None, context=None):
    """
    Fetch a variable from `context` (defaults to locals()) if defined and truthy,
    otherwise return `default`. Optionally apply `cast`.

    Args:
      name    : variable name to lookup
      default : value to return if not found or falsy
      cast    : optional function to cast the raw value
      context : dict mapping names to values (e.g. locals())

    Returns:
      The found (and cast) value or the default.
    """
    ctx = context if context is not None else locals()
    val = ctx.get(name, default)
    if val is None or val == '':
        val = default
    return cast(val) if (cast and val is not None) else val


def parse_vector(expr_str):
    """
    Parse a Feel++ vector literal into a Sympy Array.

    Accepts formats:
      '{1,2,3}', '{x,y}', or with symbol suffix ':x:y:z', e.g. '{x,y,z}:x:y:z'.

    Returns:
      Array([...]) of sympified entries.
    """
    s = str(expr_str).strip()
    m = re.match(r'^\{([^}]*)\}', s)
    if not m:
        raise ValueError(f"cannot parse vector from `{expr_str}`")
    body = m.group(1)
    parts = body.split(',') if ',' in body else body.split(':')
    return Array([sympify(p) for p in parts])

# ----------------------------------------------------------------------------
# GiNaC/C99 code emission: sympy ↔ GiNaC bridges
# ----------------------------------------------------------------------------

def syms(n):
    """Return spatial symbols [x,y,z] for dimension n."""
    return Array([x] if n == 1 else [x, y] if n == 2 else [x, y, z])


def nsyms(n):
    """Return normal symbols [nx,ny,nz] for dimension n."""
    return Array([nx] if n == 1 else [nx, ny] if n == 2 else [nx, ny, nz])


def toginac(s, symbols, rank=None):
    """
    Convert a scalar/array or Feel++ literal string into a GiNaC/C99 literal,
    appending free-symbol metadata and tensor rank.

    Args:
      s       : Sympy Expr, Array/Matrix, or string literal
      symbols : iterable of Sympy free symbols
      rank    : optional override for tensor rank (0=scalar,1=vector,2=matrix,...)

    Returns:
      String C99 literal with ':symbol...' and ':rankN' suffixes.
    """
    arr = sympify(s)
    # Handle string-encoded vectors
    if isinstance(arr, str) and arr.strip().startswith('{'):
        arr = parse_vector(arr)

    # Determine rank
    if rank is None:
        try:
            rank = Array(arr).rank()
        except Exception:
            rank = 0
    else:
        rank = int(rank)

    # Build suffixes
    symbols_list = sorted(str(sym) for sym in symbols) if symbols else []
    sym_suffix = f":{':'.join(symbols_list)}" if symbols_list else ''
    rank_suffix = f":rank{rank}"

    # Scalar literal
    if rank == 0:
        return f"{ccode(arr, standard='C99')}{sym_suffix}{rank_suffix}"

    # Array literal (vector/matrix)
    elems = arr if isinstance(arr, (Array, Matrix)) else Array(arr)
    joined = ','.join(ccode(sympify(e), standard='C99') for e in elems)
    return f"{{{joined}}}{sym_suffix}{rank_suffix}"


def sympytoginac(e):
    """
    High-level bridge: accept Expr, Array/Matrix, or literal string,
    return GiNaC/C99 literal with auto-detected rank and symbol info.
    """
    # String vector
    if isinstance(e, str) and e.strip().startswith('{'):
        vec = parse_vector(e)
        return toginac(vec, [], rank=vec.rank())

    # Sympy object
    if isinstance(e, (Expr, Array, Matrix)):
        arr = Array(e) if isinstance(e, Array) else e
        if isinstance(arr, Array) and arr.rank() == 2:
            arr = arr.tomatrix()
        sy = sympify(arr)
        rank = Array(sy).rank() if isinstance(sy, (Array, Matrix)) else 0
        syms_list = [] if not sy.free_symbols else sy.free_symbols
        return toginac(sy, syms_list, rank=rank)

    # Fallback stringify
    return toginac(str(e), [], rank=0)

# ----------------------------------------------------------------------------
# Differential operators
# ----------------------------------------------------------------------------

def dx(f):    return derive_by_array(f, [x])
def dy(f):    return derive_by_array(f, [y])
def dz(f):    return derive_by_array(f, [z])
def dt(f):    return derive_by_array(f, [t])


def grad(f, symbols):     return derive_by_array(f, symbols)
def symgrad(f, symbols):  return grad(f, symbols) + transpose(grad(f, symbols))
def div(f, symbols):      return tensorcontraction(derive_by_array(f, symbols), (0, 1))
def laplacian(f, symbols):return tensorcontraction(derive_by_array(derive_by_array(f, symbols), symbols), (0, 1))

# ----------------------------------------------------------------------------
# Tensor and flux helpers
# ----------------------------------------------------------------------------

def mult(a, b):
    a_ = Array(a)
    if a_.rank() == 0:
        return a * b
    axis = 0 if a_.rank() == 1 else 1
    return tensorcontraction(tensorproduct(a, b), (axis, axis+1))

def n(a, c=1, nsymbols=None):
    nsymbols = nsymbols or [nx, ny, nz]
    axis = 0 if Array(a).rank() == 1 else 1
    return tensorcontraction(tensorproduct(mult(c, a), nsymbols), (axis, axis+1))

def dn(a, c=1, symbols=None, nsymbols=None):
    symbols = symbols or [x, y, z]
    nsymbols = nsymbols or [nx, ny, nz]
    axis = 0 if Array(a).rank() == 0 else 1
    return tensorcontraction(tensorproduct(mult(c, grad(a, symbols)), nsymbols), (axis, axis+1))

# ----------------------------------------------------------------------------
# Additional vector/tensor operators for extended PDEs
# ----------------------------------------------------------------------------

def curl(u, symbols=[x, y, z]):
    """
    Compute the vector curl ∇×u for a 3D vector u.
    """
    return Array([
        dy(u[2]) - dz(u[1]),
        dz(u[0]) - dx(u[2]),
        dx(u[1]) - dy(u[0])
    ])


def cross(a, b):
    """
    Symbolic cross-product of two 3-vectors a and b.
    """
    return Array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ])


def hessian(f, symbols=[x, y, z]):
    """
    Return the matrix of second derivatives ∂²f/∂xi∂xj.
    """
    # diff(f, xi, xj) is a scalar Expr, so the Matrix is full of Exprs
    return Matrix([[diff(f, s1, s2) for s2 in symbols] for s1 in symbols])


def material_derivative(f, u, symbols=[x, y, z], t_symbol=t):
    """
    Material derivative D_t f = ∂_t f + u·∇f.

    Args:
      f        : scalar field (Sympy Expr)
      u        : velocity vector (Sympy Array of length len(symbols))
      symbols  : list of spatial symbols corresponding to u’s components
      t_symbol : time symbol (default `t`)

    Returns:
      Sympy Expr for D_t f.
    """
    # time derivative as a scalar Expr
    dt_f = diff(f, t_symbol)

    # spatial gradient as a list/Array of partials ∇f
    grad_f = derive_by_array(f, symbols)

    # advection term u·∇f (dot returns an Expr)
    adv = dot(u, grad_f)

    return dt_f + adv


def normal(dim):
    """
    Symbolic outward unit normal for a boundary in dimension dim.
    """
    return nsyms(dim)


def traction(sigma, dim):
    """
    Traction vector = σ · n on a boundary in dimension dim.
    """
    nvec = normal(dim)
    return mult(sigma, nvec)


def dot(a, b):
    """
    Inner product of two vectors.
    """
    return sum(a[i]*b[i] for i in range(len(a)))


def outer(a, b):
    """
    Outer product of two vectors.

    Returns
    -------
    Matrix
        A Sympy Matrix of shape (len(a), len(b)), with entry (i,j)=a[i]*b[j].
    """
    # tensorproduct gives a 2-D Array; wrapping in Matrix gives the MutableDenseMatrix tests expect
    return Matrix(tensorproduct(a, b))


def biharmonic(f, symbols=[x, y, z]):
    """
    Biharmonic operator Δ² f.
    """
    return laplacian(laplacian(f, symbols), symbols)


def integrate(f, limits):
    """
    Perform symbolic integration of `f` over given limits.

    Args:
      f      : sympy expression
      limits : list of tuples (symbol, lower, upper)

    Returns:
      Integrated expression (differed `.doit()`).
    """
    res = f
    for sym, lo, hi in limits:
        res = Integral(res, (sym, lo, hi))
    return res.doit()
