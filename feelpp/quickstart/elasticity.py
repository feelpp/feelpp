from sympy2ginac import *

#parameters={'dim':'2','mu':'1','lambda':'1','displ':{'exact','Array([x**2,0])'}};
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'Array([x,0,0])'};
#parameters={'dim':'2','mu':'1','lambda':'1','displ':'[(1/(2*pi*pi))*sin(pi*x)*cos(pi*y),(1/(2*pi*pi))*cos(pi*x)*sin(pi*y)]'}
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'[cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y) ]'}

if 'dim' in locals():
    dim=int(locals()['dim']);
else:
    dim=2

s=syms(dim);
print("s=",s);
ns=nsyms(dim);

lam1=sympify(locals()['lam1']);
print("lam1=",lam1);
lam2=sympify(locals()['lam2']);
print("Lambda=",lam2);
displ=sympify(locals()['displ']);
print("displ=",displ);

# gradient
grad_displ=grad(displ,s);
print("grad(displ)=",grad_displ);

# strain
strain =0.5 * ( grad_displ+transpose(grad_displ) );
print("strain=",strain);

# stress
I=eye(dim);
print("I=",I);
stress=Array((lam2 * tensorcontraction(strain,(0,1)) * I + 2. * lam1 *strain).tolist());
print("stress=",stress);

# surfacic forces
stressn=n(stress,1,ns);
print("stressn=",stressn);

# force density
f = -div(stress,s);
print("f=",f);

c1     = 0.5/lam1;
c2     = -lam2/(2. * lam1 * (dim*lam2 + 2.*lam1));

