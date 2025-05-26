# feelpp/sympy/codegen.py
from sympy import sympify, Array, Matrix, Expr, ccode
from .utils import parse_vector
from sympy import Matrix as SymMatrix

def toginac(s, symbols, rank=None):
    arr = sympify(s)
    if isinstance(arr, str) and arr.strip().startswith('{'):
        arr = parse_vector(arr)

    if rank is None:
        try:
            rank = Array(arr).rank()
        except:
            rank = 0
    else:
        rank = int(rank)

    symbols_list = sorted(str(sym) for sym in symbols) if symbols else []
    sym_suf = f":{':'.join(symbols_list)}" if symbols_list else ''
    rank_suf= f":rank{rank}"

    if rank == 0:
        return f"{ccode(arr,standard='C99')}{sym_suf}{rank_suf}"

    elems = arr if isinstance(arr,(Array,SymMatrix)) else Array(arr)
    joined = ','.join(ccode(e,standard='C99') for e in elems)
    return f"{{{joined}}}{sym_suf}{rank_suf}"

def sympytoginac(e):
    from sympy import Expr, Array, Matrix, sympify
    if isinstance(e,str) and e.strip().startswith('{'):
        vec = parse_vector(e)
        return toginac(vec, [], rank=vec.rank())
    if isinstance(e,(Expr,Array,Matrix)):
        arr = Array(e) if isinstance(e,Array) else e
        if isinstance(arr,Array) and arr.rank()==2:
            arr = arr.tomatrix()
        sy = sympify(arr)
        rank = Array(sy).rank() if isinstance(sy,(Array,Matrix)) else 0
        syms = [] if not sy.free_symbols else sy.free_symbols
        return toginac(sy, syms, rank=rank)
    return toginac(str(e), [], rank=0)