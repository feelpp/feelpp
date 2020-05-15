from sympy2ginac import *

if 'dim' in locals():
    dim=int(locals()['dim'])
else:
    dim=3
s=syms( dim)
ns=nsyms( dim )

if 'e' in locals():
    e=sympify(locals()['e'])
else:
    e=sympify(1)
print("e=",e)
# integrate over unit tetrahedron
if dim==3:
    I=integrate(integrate(integrate(e,(z,0,1-x-y)),(y,0,1-x)),(x,0,1))
# integrate over unit segment
if dim==1:
    I=integrate(e,(x,0,1))    
print("I=",I.evalf())
expected_value=I.evalf()





