import ukf
import numpy as np
import matplotlib.pyplot as plt

N=50

def A(x):
    return x #+ np.random.normal(0,0.1)
    
def H(x):
    return x
    
flt = ukf.Filter( dynamics = A, observe = H, defect = 10, stateDim = 1, obsDim = 1 )
stateReal = list(np.sin( np.linspace(0, 4*np.pi, N)))
signal = list(np.sin( np.linspace(0, 4*np.pi, N) ) + np.random.normal(0,0.3,N))
flt.filter(signal) 
stateEstimates = np.asmatrix(np.transpose(flt.getStateEstimateList()))
signal = np.matrix(signal)
stateReal = np.matrix(stateReal)

print("L2 norms:\n > estimation error {}\n > measurement error {}".format(np.linalg.norm(stateEstimates - stateReal),np.linalg.norm(signal - stateReal)))
#print(np.linalg.norm(stateEstimates))
#print(np.linalg.norm(stateReal))
plt.plot(np.linspace(0, 4*np.pi, N), *np.asarray(stateEstimates))
plt.plot(np.linspace(0, 4*np.pi, N), *np.asarray(signal))
# plt.plot(np.linspace(0, 4*np.pi, 100), *np.asarray(stateEstimates-stateReal))
plt.plot(np.linspace(0, 4*np.pi, N), *np.asarray(stateReal))
plt.legend(("State estimates","Observed signal","Real state"))
plt.show()
