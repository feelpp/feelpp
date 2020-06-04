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

        self.__stateEstimate = zeros(stateDim)
        self.__stateForecast = zeros(stateDim)
        self.__obsEstimate = zeros(obsDim)
        self.__obsForecast = zeros(obsDim,1+2*stateDim)
        
        self.__stateCov = zeros([stateDim,stateDim])
        self.__obsCov = zeros([obsDim,obsDim])
        self.__crossCov = zeros([stateDim,obsDim])
        self.__sigmaHk = defect*np.eye(obsDim)
        self.__gain = zeros([stateDim,obsDim])
        
        self.__sigmaScheme = defect*np.eye(stateDim)  # arbitrarily initialized with the observation defect parameter
        self.__sigmaPoints = zeros([stateDim,2*stateDim+1])
        self.__sigmaSigns = np.diag( np.concatenate( ( np.zeros(1), np.ones(stateDim), -np.ones(stateDim) ) ) )
        
        self.__weights = np.ones(2*stateDim+1)*1/6
        self.__weights[0] = 1-stateDim/3
        
    def setSigmaHk(self, M): # M is a numpy obsDim*obsDim matrix
        self.__sigmaHk = M

    def setSigmaPoints(self):
        self.__sigmaScheme = self.__sigmaSigns @ sqrtm( 3*self.__stateCov )
        self.__sigmaPoints = self.__selfEstimate * np.ones(self.stateDim) + self.__sigmaScheme        

    def getSigmaHk(self):
        return self.__sigmaHk
        
    def getSigmaPoints(self):
        return self.__sigmaPoints

    def getStateEstimate(self):
        return self.__stateEstimate

    def getStateForecast(self):
        return self.__stateForecast

    def getObsEstimate(self):
        return self.__obsEstimate

    def getObsForecast(self):
        return self.__obsForecast

    def getGain(self):
        return self.__gain

    def getWeights(self):
        return self.__weights
        
    def transformSet(self,pointSet): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.transform(pointSet[:,i])
        return pointSet

    def observeSet(self,pointSet):
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.observe(pointSet[:,i])
        return pointSet
        
    def step(self, stateEstimate, measurement):
        self.setSigmaPoints()
        self.__sigmaPoints = self.transformSet(self.__sigmaPoints)
        self.__obsForecast = self.observeSet(self.__sigmaPoints)
        
        self.__stateForecast = self.__sigmaPoints @ transpose(self.__weights)
#        self.XF = self.stateMean
        self.__stateCov = (self.__weights*(self.__sigmaPoints-self.__stateForecast)) @ transpose(self.__sigmaPoints-self.__stateForecast)
        self.__obsEstimate = self.__obsForecast @ transpose(self.__weights)
        self.__obsCov = (self.__weights*(self.__obsForecast-self.__obsEstimate)) @ transpose(self.__obsForecast-self.__obsEstimate) + self.__sigmaHk
        self.__crossCov = (self.__weights*(self.__sigmaPoints-self.__selfEstimate)) @ transpose(self.__obsForecast-self.__obsEstimate)

        self.__gain = self.__crossCov * inverse(self.__obsCov)

        self.__selfEstimate += self.__gain @ ( self.obsCurrent - self.__obsEstimate )
        
    def filter( self, measurement, maxiter = 1000, verbose = False, mode = "dynamic"):
        if mode == "dynamic":
            self.readsignal(self, measurement)
            for i in range(max(self.signal.shape)-1):
                self.X[:,i] = np.transpose(self.stateMean)
                self.forecast[:,i] = np.transpose(self.XF)
                self.step(self, mode)
                if verbose:
                    print("    sigma-points : ",self.__sigmaPoints[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.__gain)
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
                    print("    sigma-points : ",self.__sigmaPoints[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.__gain)
                    print("    last state estimate : ",self.stateMean)
                    print("    associated predicted measure : ",self.obsMean," ; real measure : ",self.signal)
                    print("    relative measure error : ",np.abs(self.obsMean-self.signal)/self.signal," ; tolerance : ",self.tol)
                i += 1
                    
            return self.stateMean
