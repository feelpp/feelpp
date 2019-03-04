import numpy as np
#from numpy import sum, mean, sqrt
#import matplotlib.pyplot as plt
from UKF import Filter
import subprocess
import time

def f(X,time,dt):
    Y = X #+ np.random.normal(0,dt)
    return Y

def h(X):
    t = time.time()
    subprocess.call("python3 /ssd/ricka/feelpp/toolboxes/pyfeelpp-toolboxes/pyfeelpptoolboxes/heat/heat.py --config-file insulating-plate.cfg --json-editions Materials.plate.k:"+str(X[0]), shell=True, stdout=-1)
    exports = open("/ssd/ricka/feelpp/feelpp/pyfeelpp/kalman/testcases/insulating-plate/np_1/np_1/heat.measures.csv","r")
    sensor = float(exports.readlines()[1])
    elapsed = time.time() - t
    print("Measurement prediction took ",elapsed," seconds")
    return sensor

setuptime = time.time()

realk = 0.2
Niter = 1000
dt = 0.1
Time = np.arange(0,Niter)
Signal = h([realk])

F = Filter
F.set(F,1,1,1/6,dt,f,h,0.1,1,tol=1e-03)
elapsed = time.time() - setuptime
print("Filter set up in ",elapsed," seconds !")

F.filter(F, Signal, maxiter = 1000, verbose=True, mode="static")

print("best estimation : ",F.Mx)

#plt.plot(Time,F.X[0,:], label="Thermal conductivity estimate")
#leg = plt.legend(loc='upper right')
#leg.get_frame().set_alpha(0.5)
#plt.show()
