import numpy as np
#from numpy import sum, mean, sqrt
#import matplotlib.pyplot as plt
#from UKF import Filter
import subprocess
from pyfeelpp import core
import sys
import pyfeelpptoolboxes.modelcore as modelcore

e=core.Environment(sys.argv,opts=modelcore.toolboxes_options("heat"))

#from pyfeelpp import discr,ts,filters
from pyfeelpptoolboxes.heat import *



def model(M,X):
    # set the parameter values
    # M.setPropertyTree({"Materials.plate.k",X[0]})
    # solve the model
    M.solve()
    # get output
    sensor=M.export("sensor")
    return sensor

f=heat(dim=2,order=1,worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()

sample = np.zeros([10,2])
sample[:,0] = np.exp(np.random.uniform(low=np.log(0.1), high=np.log(10), size=10))
sample[:,1] = np.exp(np.random.uniform(low=np.log(0.01), high=0, size=10))

for i in range(0,10):
    value = model(f,sample[i,:])
    print("*** parameter : ",sample[i,:]," temperature : ",value)
