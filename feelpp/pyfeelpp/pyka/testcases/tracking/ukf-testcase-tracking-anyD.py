import ukf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

D = 20
N = 100
T = 50

flt = ukf.Filter( dynamics = lambda x,t : x + np.random.normal(0,0.1,D)  , observe = lambda x,t : x, defect = 1e-5, stateDim = D, obsDim = D )

signal = np.cos(np.linspace(0,D*T,D*N)).reshape(D,N)

for i in range(N):
    print( "======================\nstep {}".format(i) )
    flt.step( signal[:,i] )

est = flt.getStateEstimateList()

if 1:
    #print(est)
    fig, ax = plt.subplots(D)
    for d in range(D):
        ax[d].plot(est[d,:])
        ax[d].plot(signal[d,:])
        ax[d].legend( ("component {}".format(d), "estimation {}".format(d)) )
    plt.show()

if 0:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot( xs = signal[0,:], ys = signal[1,:], zs = np.linspace(0, T, N) )
    ax.plot( xs = est[0,:], ys = est[1,:], zs = np.linspace(0, T, N) )
    plt.show()

    
