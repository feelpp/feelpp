import ukf
import advectiondiffusion as ad
import numpy as np

N = 101
dx = 0.01
dt = 0.005

def computeMatrix( mu, v ):
    return np.eye(N) - v*dt*(np.eye(N,N,1)+np.eye(N,N,-1))+mu*dt*(np.eye(N,N,1)+np.eye(N,N,-1)-2*np.eye(N))/(dx**2)

def A( x ):
    matrix = computeMatrix( x[0], x[1] )
    x[2:N+2] = matrix @ x[2:N+2]
    return x

def H( x ):
    return x[2:N+2]

flt = ukf.Filter( dynamics = A, observe = H, defect = 1, stateDim = N+2, obsDim = N )
advectionSystem = ad.Adv_Diff()

flt.step( advectionSystem.getState() )
for i in range (10):
    print(i)
    flt.step( advectionSystem.getState() )
    advectionSystem.step()

est = flt.getStateEstimateList()
for state in est:
    print(state[0:2])
