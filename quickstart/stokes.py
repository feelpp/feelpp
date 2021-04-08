from sympy2ginac import *

if 'mu' in locals():
    mu=sympify(locals()['mu']);
else:
    mu=sympify('1')
print("Mu=",mu);

if 'rho' in locals():
    rho=sympify(locals()['rho']);
else:
    rho=sympify(1)
print("rho=",rho);
nu=mu/rho
print("nu=",nu);
if 'u' in locals():
    u=sympify(locals()['u']);
else:
    u=sympify(Array([y*y-2*y+1,x*x-x]))
print("u=",u);    
if 'p' in locals():
    p=simplify(sympify(locals()['p']));
else:
    p=2*mu*(x+y-1)+1/(3*K)
print("p=",p);

if 'dim' in locals():
    dim=int(locals()['dim']);
else:
    dim=2
    
print("dim=",dim);
s=syms( dim );
ns=nsyms( dim );
div_u=div(u,s)
print("div_u=",div_u);
grad_u=grad(u,s);
print("grad_u=",grad_u);
strain=0.5 * ( grad_u+transpose(grad_u) );
print("strain=",strain);
I=eye(dim);
stress=-p * I + 2. * mu *strain;
print("stress=",stress);
stressn=n(stress,1,ns);
print("stressn=",stressn);

# compute normal component of the velocity
un=n(u,1,ns);
print("un=",un);

# compute right hand side
f = -div(stress,s);
print("f=",f);

pp = div(-2*mu*strain,s);
print("pp=",pp);
