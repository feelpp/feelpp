from sympy2ginac import *

if 'k' in locals():
    k=sympify(locals()['k']);
else:
    k=sympify(-1)
print("k=",k);
if 'p' in locals():
    p=sympify(locals()['p']);
else:
    p=x+y
print("p=",p);

if 'dim' in locals():
    dim=int(locals()['dim']);
else:
    dim=2
    
print("dim=",dim);
s=syms( dim );
ns=nsyms( dim );
grad_p=grad(p,s);
flux=-k*grad(p,s);
u=flux;
un=n(flux,1,ns);
f=div(flux,s);



