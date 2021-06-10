#from mpi4py import MPI
import numpy as np
#from numpy import sum, mean, sqrt
import matplotlib.pyplot as plt
from UKF import Filter
import subprocess
import time

#comm = MPI.COMM_WORLD
#nb_procs = comm.Get_size()
rank = 0#comm.Get_rank()

def f(X,time,dt):
    Y = X + np.random.normal(0,2*dt)
    return Y

def h(X):
    return X[0]

k = 1 
Niter = 100
dt = 0.1
Time = dt*np.arange(0,Niter)

Trajectory = k*np.cos(Time) + (1-k)*np.sin(Time)
Signal = Trajectory + np.random.normal(0, dt, size=Niter)

#comm.Bcast(Trajectory, root=0)
#comm.Bcast(Signal, root=0)

setuptime = time.time()

F = Filter
F.set(F,1,1,0.2,dt,f,h,0,1,1e-5)
print("Filter set up !")
F.filter(F, Signal, maxiter = Niter, verbose=True, mode="dynamic") 

elapsedtime = time.time() - setuptime

measurementRelativeError = np.linalg.norm(Signal[1:] - Trajectory[1:])/np.linalg.norm(Trajectory[1:])
filterRelativeError = np.linalg.norm(F.X[0] - Trajectory[1:])/np.linalg.norm(Trajectory[1:])

if (1):
    print("    ..........    Took",elapsedtime,"seconds for",F.Time,"steps")
    print(rank,"        Relative measurement error :",measurementRelativeError)
    print(rank,"        Relative filtering error :",filterRelativeError)


#print(Trajectory)
#print(Signal)
#print(F.X) 
    
#print("best estimation :",F.Mx,"performed in",elapsedtime,"seconds and",F.Time,"steps")

plt.plot(Time[:-1],F.X[0], label="estimate")
plt.plot(Time, Signal, label="signal")
plt.plot(Time, Trajectory, label="trajectory")
leg = plt.legend(loc='upper right')
leg.get_frame().set_alpha(0.5)
plt.show()
