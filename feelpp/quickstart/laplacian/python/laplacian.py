from sympy2ginac import *

if 'compute_pde_coefficients' in locals() and locals()['compute_pde_coefficients']:
    compute_pde_coefficients=locals()['compute_pde_coefficients'];
else:
    compute_pde_coefficients='true'

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
J=tensorcontraction(tensorproduct(u,u),(0,1))/k;

if compute_pde_coefficients=='true':
    un=n(flux,1,ns);
    f=div(flux,s);
    g=p
    r_2=un-r_1*p
else:
    if 'un' in locals() and locals()['un']:
        un=sympify(locals()['un'])
    else:
        raise RuntimeError('invalid data, un was not specified')    
    if 'f' in locals() and locals()['f']:
        f=sympify(locals()['f'])
    else:
        raise RuntimeError('invalid data, f was not specified')         
    if 'g' in locals() and locals()['g']:
        g=sympify(locals()['g'])
    else:
        raise RuntimeError('invalid data, g was not specified')                 
    if 'r_2' in locals() and locals()['r_2']:
        r_2=sympify(locals()['r_2'])
    else:
        raise RuntimeError('invalid data, r_2 was not specified')        
