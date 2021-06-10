import ukf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 100
T = 20

flt = ukf.Filter( dynamics = lambda x : x + np.random.normal(0,0.1,2)  , observe = lambda x : x, defect = 1e-5, stateDim = 2, obsDim = 2 )

signal = np.cos(np.linspace(0,2*T,2*N)).reshape(2,N)

for i in range(N):
    flt.step( signal[:,i] )

est = flt.getStateEstimateList()

if 1:
    #print(est)
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot( signal[0,:] )
    ax1.plot( est[0,:] )
    ax1.legend( ("1st component", "1st component estimation") )
    ax2.plot( signal[1,:] )
    ax2.plot( est[1,:] )
    ax2.legend( ("2nd component", "2nd component estimation") )
    plt.show()

if 0:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot( xs = signal[0,:], ys = signal[1,:], zs = np.linspace(0, T, N) )
    ax.plot( xs = est[0,:], ys = est[1,:], zs = np.linspace(0, T, N) )
    plt.show()

    
