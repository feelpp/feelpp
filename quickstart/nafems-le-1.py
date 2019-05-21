from sympy2ginac import *
from sympy.physics.units import *

#parameters={'dim':'2','mu':'1','lambda':'1','displ':{'exact','Array([x**2,0])'}};
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'Array([x,0,0])'};
#parameters={'dim':'2','mu':'1','lambda':'1','displ':'[(1/(2*pi*pi))*sin(pi*x)*cos(pi*y),(1/(2*pi*pi))*cos(pi*x)*sin(pi*y)]'}
#parameters={'dim':'3','mu':'1','lambda':'1','displ':'[cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y) ]'}

dim=2;
s=syms(dim);
ns=nsyms(dim);


E=210*(10**9);
nu=0.3;
rho=7800;
stressn=dict();
stressn['BC']=10*(10**6);

Lambda=(E*nu)/((1+nu)*(1-2*nu));
Mu=E/(2*(1+nu));

displ_x=dict();displ_y=dict();
displ_y['DC']=sympify(0)
displ_x['AB']=sympify(0)
print("displ_x=",displ_x);
print("displ_y=",displ_y);


# force density
f=dict();
f['material'] = sympify(0);
print("f=",f);

c1=dict();c2=dict();
c1['material']     = 0.5/Mu;
c2['material']     = -Lambda/(2. * Mu * (dim*Lambda + 2.*Mu));
print("c1=",c1," c2=", c2);
