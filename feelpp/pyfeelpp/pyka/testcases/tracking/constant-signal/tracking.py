import numpy as np

#import matplotlib.pyplot as plt
#from matplotlib2tikz import save as tikz_save

from base import *
from UKF import Filter
import subprocess
import time

def f(X,time,dt):
    Y = X + np.random.normal(0,1)
    return Y

def h(X):
    return X[0]

realk = 10
Niter = 10000
dt = 1
Time = np.arange(0,Niter)
#Signal = h([realk])*np.ones(Niter)

setuptime = time.time()

F = Filter
F.set(F,1,1,1/2,dt,f,h,0.1,1,1e-5)
#F.readsignal(F,Signal)
print("Filter set up !")
F.filter(F, realk, maxiter = Niter, verbose=False, mode="static") 

elapsedtime = time.time() - setuptime

print("best estimation :",F.Mx,"performed in",elapsedtime,"seconds and",F.Time,"steps")

#plt.plot(Time,F.X[0,:], label="Thermal conductivity estimate")
#leg = plt.legend(loc='upper right')
#leg.get_frame().set_alpha(0.5)
#plt.show()
#tikz_save("fig.tex")

export = open("export","w")
for i in range(F.Time):
    if i%60==0:
        export.write("("+str(Time[i])+", "+str(F.X[0,i])+")\n")
export.close()

