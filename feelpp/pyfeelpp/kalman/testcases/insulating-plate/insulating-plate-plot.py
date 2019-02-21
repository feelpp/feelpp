import numpy as np
#from numpy import sum, mean, sqrt
#import matplotlib.pyplot as plt
#from UKF import Filter
import subprocess

def model(X):
#    print("current state to observe : "+str(X))
    subprocess.call("python3 /ssd/ricka/feelpp/toolboxes/pyfeelpp-toolboxes/pyfeelpptoolboxes/heat/heat.py --config-file insulating-plate.cfg --json-editions Materials.plate.k:"+str(X[0])+" BoundaryConditions.temperature.Robin.gamma2.expr1:"+str(X[1]), shell=True)
    exports = open("/ssd/ricka/feelpp/feelpp/pyfeelpp/kalman/testcases/insulating-plate/np_1/np_1/heat.measures.csv","r")
    sensor = float(exports.readlines()[1])
    print("current observation : "+str(sensor))
    return sensor

sample = np.zeros([10,2])
sample[:,0] = np.exp(np.random.uniform(low=np.log(0.1), high=np.log(10), size=10))
sample[:,1] = np.exp(np.random.uniform(low=np.log(0.01), high=0, size=10))

for i in range(0,10):
    value = model(sample[i,:])
    print("*** parameter : ",sample[i,:]," temperature : ",value)
