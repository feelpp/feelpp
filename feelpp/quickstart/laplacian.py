from sympy2ginac import *

if 'k' in locals():
    k=sympify(locals()['k']);
else:
    k=sympify(1)
print("k=",k);
if 'r_1' in locals():
    r_1=sympify(locals()['r_1']);
else:
    r_1=sympify(1)

print("r_1=",r_1);
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
if 'r_2' in locals() and locals()['r_2']:
    r_2=sympify(locals()['r_2']);
else:
    # right hand side Robin condition when given r_1
    # r_1 is positive or zero and the flux un=r_1*p+r_2
    r_2=un-r_1*p
print("r_2=",r_2)


