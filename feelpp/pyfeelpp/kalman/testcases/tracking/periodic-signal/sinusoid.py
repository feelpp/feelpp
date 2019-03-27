#from mpi4py import MPI
import numpy as np

import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

from UKF import Filter
import subprocess
import time

#comm = MPI.COMM_WORLD
#nb_procs = comm.Get_size()
rank = 0#comm.Get_rank()

# FUNCTION DEFINITIONS

def f(X,time,dt):
    Y = X + 0.5*np.sin(time*dt) + 0.5*np.cos(time*dt)
    return Y

def h(X):
    return X[0]

# SETUP

k = 0.6
Niter = 100
dt = 0.1
Time = dt*np.arange(0,Niter)

Trajectory = np.zeros(Niter)
for i in range(Niter-1):
    ki = np.random.normal(k,0.1)
    Trajectory[i+1] = Trajectory[i] + ki*np.sin(Time[i]) + (1-ki)*np.cos(Time[i])
Signal = Trajectory + np.random.normal(0, 1, size=Niter)

#comm.Bcast(Trajectory, root=0)
#comm.Bcast(Signal, root=0)

setuptime = time.time()

F = Filter
F.set(F,1,1,0,dt,f,h,0.1,1,1e-5)
print("Filter set up !")
F.filter(F, Signal, maxiter = Niter, verbose=False, mode="dynamic") 

elapsedtime = time.time() - setuptime

measurementRelativeError = np.linalg.norm(Signal[1:] - Trajectory[1:])/np.linalg.norm(Trajectory[1:])
filterRelativeError = np.linalg.norm(F.X[0] - Trajectory[1:])/np.linalg.norm(Trajectory[1:])

if (1):
    print("    ..........    Took",elapsedtime,"seconds for",F.Time,"steps")
    print("        Relative measurement error :",measurementRelativeError)
    print("        Relative filtering error :",filterRelativeError)


#print(Trajectory)
#print(Signal)
#print(F.X) 

plt.plot(Time[0:-1],F.X[0], label="estimate")
plt.plot(Time, Signal, label="signal")
plt.plot(Time, Trajectory, label="trajectory")
leg = plt.legend(loc='upper right')
leg.get_frame().set_alpha(0.5)
plt.show()
tikz_save("sinusoid.tex")
