from sympy2ginac import *

#parameters={'dim':'2','mu':'1','lambda':'1','displ':{'exact','Array([x**2,0])'}};
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'Array([x,0,0])'};
#parameters={'dim':'2','mu':'1','lambda':'1','displ':'[(1/(2*pi*pi))*sin(pi*x)*cos(pi*y),(1/(2*pi*pi))*cos(pi*x)*sin(pi*y)]'}
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'[cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y) ]'}

dim=int(parameters['dim']);
s=syms(dim);
print("s=",s);
ns=nsyms(dim);

Mu=sympify(parameters['mu']);
Lambda=sympify(parameters['lambda']);

displ=dict();
displ['exact']=sympify(Array([x**2,0]));
displ['Dirichlet']=displ['exact'];
print("displ=",displ);

# gradient
grad_displ=dict();
grad_displ['exact']=grad(displ['exact'],s);
print("grad(displ)=",grad_displ);

# strain
strain=dict();
strain['exact'] =0.5 * ( grad_displ['exact']+transpose(grad_displ['exact']) );
print("strain=",strain);

# stress
stress=dict();
I=eye(dim);
stress['exact']=Lambda * tensorcontraction(strain['exact'],(0,1)) * I + 2. * Mu *strain['exact'];
print("stress=",stress);

# surfacic forces
stressn=dict();
stressn['Neumann']=n(stress['exact'],1,ns);
print("stressn=",stressn);

# force density
f=dict();
f['exact'] = -div(stress['exact'],s);
print("f=",f);

c1['exact']     = 0.5/Mu;
c2['exact']     = -Lambda/(2. * Mu * (dim*Lambda + 2.*Mu));

