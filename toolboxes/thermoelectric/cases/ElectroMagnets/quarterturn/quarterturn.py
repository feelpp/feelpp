#------------------------------------------------------------------------------
# Compute the exact solution of Heat Equation on Torus
#------------------------------------------------------------------------------

from matplotlib import cm
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.linalg import solve
import pint
import pint_pandas
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.setup_matplotlib(True)
Q_ = ureg.Quantity
#-------------------------#------------------------------------------------------------------------------
# Compute the exact solution of Heat Equation on Torus
#------------------------------------------------------------------------------


#--------------------------------------
# Parameters

#%% units

def computeT(r,scaledlength,length_unit):    
    rint = 75*scaledlength     # m
    rext = 100.2*scaledlength  # m
    z1 = 50/2.*scaledlength    # m

    h = 80000*ureg.W/(ureg.m**2*ureg.K)       # W/m2/K
    T_c = 293*ureg.K       # K
    U = 1. * ureg.V        # V
    k = 380*ureg.W/(ureg.K*ureg.m)          # W/m/K
    sigma = 58.e+6*ureg.S/ureg.m   # S.m-1
    sigma.to_base_units()
    n = 1000


    #--------------------------------------
    # Compute of constants

    r1 = rint
    r2 = rext
    h1 = h
    h2 = h
    Tw1 = T_c
    Tw2 = T_c

    sigma0 = sigma
    V0 = U

    a = sigma0/(2*k)*(V0/(2*math.pi))**2
    a.check('[temperature]')
    b = k*(1/(h1*r1)+1/(h2*r2))+math.log(r2/r1)
    b.check('') # dimensionless
    c = math.log(r2/r1)*math.log(r2.magnitude*r1.magnitude)+2*k * \
        (math.log(r1.magnitude)/(h1*r1)+math.log(r2.magnitude)/(h2*r2))
    d=((Tw2-Tw1)/b+a*c/b)/(2*a)
    d.check('[length]')
    r0 = math.exp(d.magnitude)*length_unit
    r0.check('[length]')
    Tm = 2*a*k/(h1*r1+h2*r2)*math.log(r2/r1)
    Tm += (h1*r1*Tw1+h2*r2*Tw2)/(h1*r1+h2*r2)
    Tm += a*(h1*r1*math.log(r1/r0)**2+h2*r2*math.log(r2/r0)**2)/(h1*r1+h2*r2)
    Tm=Tm.to_base_units()
    Tm.check('[temperature]')
    print("a = {}".format(a))
    print("b = {}".format(b))
    print("c = {}".format(c))
    print("r0 = {}".format(r0))
    print("r1 = {}".format(r1))
    print("r2 = {}".format(r2))
    print("Tmax = {}".format(Tm))

    I = sigma * U/(2*math.pi) * math.log(r2/r1) * 2*z1
    I=I.to_base_units()
    I.check('[current]')
    print("I = {0}".format(I))
    T=-a*np.log(r/r0)**2+Tm
    T.check('[temperature]')
    return T

#LengthUnit = ureg.m
#ScaledUnit = 1.e-3*ureg.m
LengthUnit = ureg.mm
ScaledUnit = 1*ureg.mm

fig = plt.figure()
ax = plt.gca()

r = np.arange(start=r1.magnitude, stop=r2.magnitude, step=(r2-r1).magnitude/10)*LengthUnit
#print(r)
T = computeT(r, ScaledUnit, LengthUnit).to_base_units()
print(T)
plt.plot(r, T, label='Tm', marker='x')

ax.set_title('Maximum temperature vs radius')
plt.show()

