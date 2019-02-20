import numpy as np
#from numpy import sum, mean, sqrt
#import matplotlib.pyplot as plt
from UKF import Filter
import subprocess

def f(X,time,dt):
    Y = X + np.random.normal(0,dt)
    return Y

def h(numpyarray):
    sensor = np.zeros(len(numpyarray))
    for k in range(len(numpyarray)):
        subprocess.call("python3 /ssd/ricka/feelpp/toolboxes/pyfeelpp-toolboxes/pyfeelpptoolboxes/heat/heat.py --config-file insulating-plate.cfg --json-editions Materials.plate.k:"+str(numpyarray[k]), shell=True)
        exports = open("/ssd/ricka/feelpp/feelpp/pyfeelpp/kalman/testcases/insulating-plate/np_1/np_1/heat.measures.csv","r")
        sensor[k] = float(exports.readlines()[1])
    return sensor

realk = 0.2
Niter = 1000
dt = 0.01
Time = np.arange(0,Niter)
Signal = h([0.2])*np.ones(Niter)

F = Filter
F.set(F,1,1,0.2,dt,f,h,0.01)
F.readsignal(F,Signal)
print("Filter set up !")
F.filter(F)

print(F.X[0,-1])

#plt.plot(Time,F.X[0,:], label="Thermal conductivity estimate")
#leg = plt.legend(loc='upper right')
#leg.get_frame().set_alpha(0.5)
#plt.show()
