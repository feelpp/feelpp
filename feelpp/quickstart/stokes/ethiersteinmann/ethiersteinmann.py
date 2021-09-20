import sympy2ginac
from sympy2ginac import *
from sympy import *
import  math
#parameters={'dim':'2','mu':'1','lambda':'1','velocity':{'exact','Array([x**2,0])'}};
#parameters={'dim':'3','mu':'1','lambda':'1','velocity':'Array([x,0,0])'};
#parameters={'dim':'2','mu':'1','lambda':'1','velocity':'[(1/(2*pi*pi))*sin(pi*x)*cos(pi*y),(1/(2*pi*pi))*cos(pi*x)*sin(pi*y)]'}
#parameters={'dim':'3','mu':'1','lambda':'1','velocity':'[cos(Pi*x)*cos(Pi*y)*cos(Pi*z), cos(Pi*y)*sin(Pi*x)*sin(Pi*z), cos(Pi*x)*cos(Pi*z)*sin(Pi*y) ]'}

if 'dim' in locals():
    dim=int(locals()['dim']);
else:
    dim=2

if 'exact' in locals():
    exact=True if locals()['exact']=="1" else False;
else:
    exact=True

s=syms(dim);
print("s=",s);
ns=nsyms(dim);

if 'mu' in locals():
    mu=sympify(locals()['mu']).evalf()
else:
    mu=0.1    
print("mu=",mu);

potential = - (pi/4)*(pi/4)*0.5 * ( exp(2*(pi/4)*x) + exp(2*(pi/4)*y) + exp(2*(pi/4)*z) + 2*sin((pi/4)*x+(pi/2)*y)*cos((pi/4)*z+(pi/2)*x)*exp((pi/4)*(y+z)) + 2*sin((pi/4)*y+(pi/2)*z)*cos((pi/4)*x+(pi/2)*y)*exp((pi/4)*(z+x)) + 2*sin((pi/4)*z+(pi/2)*x)*cos((pi/4)*y+(pi/2)*z)*exp((pi/4)*(x+y)))

print("-- p=",potential)
velocity = Array([-(pi/4)*(exp((pi/4)*x)*sin((pi/4)*y+(pi/2)*z)+exp((pi/4)*z)*cos((pi/4)*x+(pi/2)*y)),
                  -(pi/4)*(exp((pi/4)*y)*sin((pi/4)*z+(pi/2)*x)+exp((pi/4)*x)*cos((pi/4)*y+(pi/2)*z)),
                  -(pi/4)*(exp((pi/4)*z)*sin((pi/4)*x+(pi/2)*y)+exp((pi/4)*y)*cos((pi/4)*z+(pi/2)*x))])
print("-- u=",velocity)

if exact:
    # gradient
    grad_velocity=grad(velocity,s)
    print("grad(velocity)=",grad_velocity);
    
    # strain
    strain = 0.5 * ( grad_velocity+transpose(grad_velocity) )
    print("strain=",strain);
    
    # stress
    I=eye(dim);
    print("I=",I);
    stress=- potential * I + 2. * mu *strain;
    print("stress=",stress);
    
    # surfacic forces
    stressn=n(stress,1,ns)
    
    # force density
    f = -sympy2ginac.div(stress,s)
else:
    if 'stressn' in locals():
        stressn=sympify(locals()['stressn'])
    else:
        stressn=sympify('0')
    f=sympify(locals()['f']);

print("stressn=",stressn);
print("f=",f);


