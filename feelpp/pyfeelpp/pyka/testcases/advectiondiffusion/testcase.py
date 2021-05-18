import ukf
import advectiondiffusion as ad
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

mu, v = 0.05, 0.5
a, b = 1, 0.25
omega = 10
N = 101
dx = 1/(N-1)
dt = 0.15 # min( dx**2/(2*mu), dx/mu )/50

def computeMatrix( v, mu ):
    M = np.eye(N) - v*dt*(np.eye(N,N,1)+np.eye(N,N,-1))+mu*dt*(np.eye(N,N,1)+np.eye(N,N,-1)-2*np.eye(N))/(dx**2)
    return M

def computeCovariance( theta ):
    return np.diag(np.concatenate((0.125*theta,np.zeros(N)))) # 12.5% of given initial values, model assumed exact

def A( x, t ):
    x[5:] = computeMatrix( x[0], x[1] ) @ x[5:]
    x[5] = x[2] + x[3]*np.sin( x[4]*t )
    x[-1] = 2*x[-2] - x[-3]
    x[0:5] += np.random.normal( 0, 1e-8, 5 ) # random walk for constant unknown parameters
    return x

def H( x, t ):
    return x[[9,29,49,84]]

def h( x, t ):
    return x[[5,25,45,80]]

initialState = np.ones((N+5,1))
#initialState[0:4,0] = [1,1,a,b]
flt = ukf.Filter( dynamics = A, observe = H, defect = 0.002, stateDim = N+5, obsDim = 4 )
flt.setState( initialState )
#flt.setSigmaHk( np.zeros((N,N)) )
advectionSystem = ad.Adv_Diff( N = N, dx = dx, dt = dt, v=v, mu = mu, a = a, b = b )
advectionSystem._state = np.ones(N) # as prescribed in Lal

flt.step( h(advectionSystem.getState(), 0) )
for i in range (667):
    print(i)
    print("observation : {}".format(h(advectionSystem.getState(), flt.time)))
    flt.step( h(advectionSystem.getState(), flt.time) )
    advectionSystem.step()
    advectionSystem.setBC( a + b*np.sin( omega*i*dt ), 2*advectionSystem.getState()[-2] - advectionSystem.getState()[-3] )
    advectionSystem.resetBC()
    print( flt.getStateEstimate()[0:5] )
    
est = flt.getStateEstimateList()
#for state in est:
#    print(state[0:2])
print(np.linalg.norm( est[2:,-1]-advectionSystem.getState() ))
plt.plot( est[2:,-1] )
plt.plot( advectionSystem.getState() )
plt.show()
