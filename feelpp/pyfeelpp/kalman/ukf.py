import numpy as np
from numpy import mean, sum, dot, linalg, zeros, ones, transpose
from numpy.linalg import inv, cholesky
from scipy.linalg import sqrtm
from base import *

class Filter:

    def __init__(self, dynamics, observe, defect, stateDim, obsDim):
        self.observe = observation
        self.transform = dynamics
        self.stateDim = stateDim
        self.obsDim = obsDim

        self.stateEstimate = zeros(stateDim)
        self.stateForecast = zeros(stateDim)
        self.obsEstimate = zeros(obsDim)
        self.obsForecast = zeros(obsDim,1+2*stateDim)
        
        self.stateCov = zeros([stateDim,stateDim])
        self.obsCov = zeros([obsDim,obsDim])
        self.crossCov = zeros([stateDim,obsDim])
        self.sigmaHk = defect*np.eye(obsDim)
        self.gain = zeros([stateDim,obsDim])
        
        self.sigmaScheme = defect*np.eye(stateDim)  # arbitrarily initialized with the observation defect parameter
        self.sigmaPoints = zeros([stateDim,2*stateDim+1])
        self.sigmaSigns = np.diag( np.concatenate( ( np.zeros(1), np.ones(stateDim), -np.ones(stateDim) ) ) )
        
        self.weights = np.ones(2*stateDim+1)*1/6
        self.weights[0] = 1-stateDim/3

    def transformSet(self,pointSet): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.transform(pointSet[:,i])
        return pointSet

    def observeSet(self,pointSet):
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.observe(pointSet[:,i])
        return pointSet
        
    def setSigmaHk(self, M): # M is a numpy obsDim*obsDim matrix
        self.sigmaHk = M

    def setSigmaPoints(self):
        self.sigmaScheme = self.sigmaSigns @ sqrtm( 3*self.stateCov )
        self.sigmaPoints = self.stateEstimate * np.ones(self.stateDim) + self.sigmaScheme
        
    def step(self, stateEstimate, measurement):
        self.setSigmaPoints()
        self.sigmaPoints = self.transformSet(self.sigmaPoints)
        self.obsForecast = self.observeSet(self.sigmaPoints)
        
        self.stateForecast = self.sigmaPoints @ transpose(self.weights)
#        self.XF = self.stateMean
        self.stateCov = (self.weights*(self.sigmaPoints-self.stateForecast)) @ transpose(self.sigmaPoints-self.stateForecast)
        self.obsEstimate = self.obsForecast @ transpose(self.weights)
        self.obsCov = (self.weights*(self.obsForecast-self.obsEstimate)) @ transpose(self.obsForecast-self.obsEstimate) + self.sigmaHk
        self.crossCov = (self.weights*(self.sigmaPoints-self.stateEstimate)) @ transpose(self.obsForecast-self.obsEstimate)

        self.gain = self.crossCov * inverse(self.obsCov)

        self.stateEstimate += self.gain @ ( self.obsCurrent - self.obsEstimate )
        
    def filter( self, measurement, maxiter = 1000, verbose = False, mode = "dynamic"):
        if mode == "dynamic":
            self.readsignal(self, measurement)
            for i in range(max(self.signal.shape)-1):
                self.X[:,i] = np.transpose(self.stateMean)
                self.forecast[:,i] = np.transpose(self.XF)
                self.step(self, mode)
                if verbose:
                    print("    sigma-points : ",self.sigmaPoints[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.gain)
                    print("    last state estimate : ",self.stateMean)
                    print("    associated predicted measure : ",self.obsMean," ; real measure : ",self.signal[i])
                    print("    relative measure error : ",np.abs(self.obsMean-self.signal[i])/self.signal[i]," ; tolerance : ",self.tol)
#                    if np.abs(self.obsMean-self.signal[i])/self.signal[i] < self.tol or i > maxiter:
#                        return self.stateMean
                    
        elif mode == "static":
            i = 0   
            self.signal = measurement
            self.X = zeros([self.stateDim,maxiter]) # KEEPS TRACK OF THE STATES BEST ESTIMATE
            self.forecast = zeros([self.stateDim,maxiter]) # KEEPS TRACK OF THE FORECAST
            
            while np.abs(self.obsMean-self.signal)/self.signal > self.tol and i < maxiter:
                self.X[:,i] = np.transpose(self.stateMean)
                self.forecast[:,i] = np.transpose(self.XF)
                self.step(self, mode)
                if verbose:
                    print("    sigma-points : ",self.sigmaPoints[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.gain)
                    print("    last state estimate : ",self.stateMean)
                    print("    associated predicted measure : ",self.obsMean," ; real measure : ",self.signal)
                    print("    relative measure error : ",np.abs(self.obsMean-self.signal)/self.signal," ; tolerance : ",self.tol)
                i += 1
                    
            return self.stateMean
