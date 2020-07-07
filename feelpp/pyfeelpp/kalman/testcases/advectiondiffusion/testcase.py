import ukf
import advectiondiffusion as ad
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

N = 101
dx = 1/(N-1)
dt = 0.000000001

def computeMatrix( mu, v ):
    M = np.eye(N) - v*dt*(np.eye(N,N,1)+np.eye(N,N,-1))+mu*dt*(np.eye(N,N,1)+np.eye(N,N,-1)-2*np.eye(N))/(dx**2)
    return M

def A( x ):
    BC = [2,1] # x[[2,-1]]
    matrix = computeMatrix( x[0], x[1] )
    x[2:] = matrix @ x[2:]
    x[[2,-1]] = BC
    return x

def H( x ):
    return x[2:N+2]

initialState = np.zeros((N+2,1))
initialState[[0,1,2,-1],0] = [1,1,2,1]
flt = ukf.Filter( dynamics = A, observe = H, defect = 0.0001, stateDim = N+2, obsDim = N )
flt.setState( initialState )
#flt.setSigmaHk( np.zeros((N,N)) )
advectionSystem = ad.Adv_Diff( N = N, dx = dx, dt = dt, v=1, mu = 0.01, a=1, b=0 )

flt.step( advectionSystem.getState() )
for i in range (50):
    print(i)
    flt.step( advectionSystem.getState() )
    advectionSystem.step()
    print( flt.getStateEstimate()[0:2] )
    
est = flt.getStateEstimateList()
#for state in est:
#    print(state[0:2])
print(np.linalg.norm( est[2:,-1]-advectionSystem.getState() ))
plt.plot( est[2:,-1] )
plt.plot( advectionSystem.getState() )
plt.show()
