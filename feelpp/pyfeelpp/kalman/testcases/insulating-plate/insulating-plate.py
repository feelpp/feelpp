import numpy as np
#from numpy import sum, mean, sqrt
#import matplotlib.pyplot as plt
from UKF import Filter
import subprocess

def f(X,time,dt):
    Y = X + np.random.normal(0,dt)
#    print("sigma-points "+str(F.SigPts))
    print("estimation "+str(np.mean(F.SigPts)))
    return Y

def h(X):
#    print("current state to observe : "+str(X))
    subprocess.call("python3 /ssd/ricka/feelpp/toolboxes/pyfeelpp-toolboxes/pyfeelpptoolboxes/heat/heat.py --config-file insulating-plate.cfg --json-editions Materials.plate.k:"+str(X[0]), shell=True)
    exports = open("/ssd/ricka/feelpp/feelpp/pyfeelpp/kalman/testcases/insulating-plate/np_1/np_1/heat.measures.csv","r")
    sensor = float(exports.readlines()[1])
    print("current observation : "+str(sensor))
    return sensor

realk = 0.2
Niter = 1000
dt = 0.1
Time = np.arange(0,Niter)
Signal = h([realk])*np.ones(Niter)

F = Filter
F.set(F,1,1,2/3,dt,f,h,0.1,0.21)
F.readsignal(F,Signal)
print("Filter set up !")
F.filter(F)

print("best estimation : ",F.Mx)

#plt.plot(Time,F.X[0,:], label="Thermal conductivity estimate")
#leg = plt.legend(loc='upper right')
#leg.get_frame().set_alpha(0.5)
#plt.show()
