from sympy import *
x, y, z, t = symbols("x y z t")
nx, ny, nz = symbols("nx ny nz")
alpha, T = symbols("alpha T")

def syms(n):
    return Array( [x] if n==1 else [x,y] if n==2 else [x,y,z] );

def nsyms(n):
    return Array( [nx] if n==1 else [nx,ny] if n==2 else [nx,ny,nz] )

def toginac(s,symbols):
    #print("s=",sympify(s));
    strsymbols=[]
    for i in symbols:
        strsymbols.append(str(i))
    strsymbols=sorted(strsymbols);
    if Array(sympify(s)).rank() == 0:
        if strsymbols:
            return str(ccode(sympify(s),standard='C99'))+':' + ':'.join(str(e) for e in strsymbols);
        else:
            return str(ccode(sympify(s),standard='C99'))
    if strsymbols:
        return '{' + ','.join(str(ccode(sympify(e),standard='C99')) for e in s) + '}:' + ':'.join(str(e) for e in strsymbols)
    else:
        return '{' + ','.join(str(ccode(sympify(e),standard='C99')) for e in s) + '}'
        

def sympytoginac(e):
    if isinstance(e,Expr) or isinstance(e,Array):
        return toginac(sympify( e ), [] if len( e.free_symbols)==0 else e.free_symbols );
    return str(e);

def dx(f):
    return derive_by_array(f,[x]);

def dy(f):
    return derive_by_array(f,[y]);

def dz(f):
    return derive_by_array(f,[z]);

def dt(f):
    return derive_by_array(f,[t]);

def grad(f,symbols):
    return derive_by_array(f,symbols);

def symgrad(f,symbols):
    a=grad(f,symbols);
    return a+transpose(a);

def div(f,symbols):
    return tensorcontraction(derive_by_array(f,symbols), (0,1))

def laplacian(f,symbols):
    return tensorcontraction(derive_by_array(derive_by_array(f,symbols),symbols), (0,1))

def mult(a,b):
    a_=Array(a);
    if a_.rank()==0 :
        return a*b;
    else:
        axe=0 if Array(a).rank()==1 else 1;
        return tensorcontraction(tensorproduct(a,b),(axe,axe+1));

def n(a,c=1,nsymbols=[nx,ny,nz]):
    axe=0 if Array(a).rank()==1 else 1;
    #print("axe:",axe);
    #print("a:",a);
    return tensorcontraction(tensorproduct(mult(c,a),nsymbols),(axe,axe+1));

def dn(a,c=1,symbols=[x,y,z],nsymbols=[nx,ny,nz]):
    axe=0 if Array(a).rank()==0 else 1;
    return tensorcontraction(tensorproduct(mult(c,grad(a,symbols)),nsymbols),(axe,axe+1));
