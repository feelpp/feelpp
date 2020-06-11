import ukf
import numpy as np

def f(x):
    return x + np.random.normal(0,0.1)
    
def h(x):
    return x
    
F = ukf.Filter(f, h, 1000, 1, 1)
signal = list(np.sin( np.linspace(0, 4*np.pi, 100) ))
F.filter(signal) 
stateEstimates = np.asmatrix(np.transpose(F.getStateEstimateList()))
signal = np.matrix(signal)

print(np.linalg.norm(stateEstimates - signal))
