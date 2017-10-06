from sympy2ginac import *
init_printing()
#parameters={'dim':'2','mu':'1','lambda':'1','displ':'[x*y,y^2]'};
parameters={'dim':'2','mu':'1','lambda':'1','displ':'[(1/(2*pi*pi))*sin(pi*x)*cos(pi*y),(1/(2*pi*pi))*cos(pi*x)*sin(pi*y)]'}
parameters={'dim':'3','mu':'1','lambda':'1','displ':'[cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y) ]'}

dim=int(parameters['dim']);
s=syms(dim);
print("s=",s);
ns=nsyms(dim);

Mu=sympify(parameters['mu']);
Lambda=sympify(parameters['lambda']);

displ=sympify(parameters['displ']);
print("displ=",displ);
print("grad(displ)=",grad(displ,s));
strain=0.5 * ( grad(displ,s)+transpose(grad(displ,s)) );
print("strain=",strain);
I=eye(dim);
print("I=",I);
stress=Lambda * tensorcontraction(strain,(0,1)) * I + 2. * Mu *strain;
print("stress=",stress);
f = div(stress,s);
print("f=",f);
c1     = 0.5/Mu;
c2     = -Lambda/(2. * Mu * (dim*Lambda + 2.*Mu));

