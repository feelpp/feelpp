from sympy2ginac import *

if 'k' in locals() and locals()['k']:
    k=sympify(locals()['k']);
else:
    k=sympify(1)

if 'r_1' in locals() and locals()['r_1']:
    r_1=sympify(locals()['r_1']);
else:
    r_1=sympify(1)


if 'p' in locals():
    p=sympify(locals()['p']);
else:
    p=x+y

if 'dim' in locals():
    dim=int(locals()['dim']);
else:
    dim=2
    

s=syms( dim );
ns=nsyms( dim );
grad_p=grad(p,s);
flux=-k*grad(p,s);
u=flux;
if 'un' in locals() and locals()['un']:
    un=sympify(locals()['un'])
else:
    un=n(flux,1,ns);
if 'f' in locals() and locals()['f']:
    f=sympify(locals()['f'])
else:
    f=div(flux,s);
if 'g' in locals() and locals()['g']:
    g=sympify(locals()['g'])
else:
    g=p
if 'r_2' in locals() and locals()['r_2']:
    r_2=sympify(locals()['r_2']);
else:
    # right hand side Robin condition when given r_1
    # r_1 is positive or zero and the flux un=r_1*p+r_2
    r_2=un-r_1*p



