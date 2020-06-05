import numpy as np
from numpy import mean, sum, dot, linalg, zeros, ones, transpose
from numpy.linalg import inv, cholesky
from scipy.linalg import sqrtm
from base import *

class Filter:

    def __init__( self, dynamics, observe, defect, stateDim, obsDim ):
        self.observe = observation
        self.transform = dynamics
        
        self.stateDim = stateDim
        self.obsDim = obsDim

        self.__stateEstimateList = []
        self.__obsEstimateList = []
        
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
        
    def setSigmaHk( self, M ): # M is a numpy obsDim*obsDim matrix
        self.__sigmaHk = M

    def setSigmaPoints( self ):
        self.__sigmaScheme = self.__sigmaSigns @ sqrtm( 3*self.__stateCov )
        self.__sigmaPoints = self.__stateEstimate * np.ones((1,self.stateDim)) + self.__sigmaScheme        

    def __setStateEstimate( self, value ):
        self.__stateEstimate = value
        
    def getSigmaHk( self ):
        return self.__sigmaHk
        
    def getSigmaPoints( self ):
        return self.__sigmaPoints

    def getStateEstimate( self ):
        return self.__stateEstimate

    def getStateForecast( self ):
        return self.__stateForecast

    def getObsEstimate( self ): # this is the weighted mean of H(A(sigma-points))
        return self.__obsEstimate

    def getObsForecast( self ): # this is the set of H(A(sigma-points))
        return self.__obsForecast

    def getGain( self ):
        return self.__gain

    def getWeights( self ):
        return self.__weights

    def getStateEstimateList( self ):
        return self.__stateEstimateList

    def getObsEstimateList( self ):
        return self.__obsEstimateList

    def getNiter( self ):
        return len(self.getStateEstimateList())
        
    def transformSet( self, pointSet ): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.transform(pointSet[:,i])
        return pointSet

    def observeSet( self, pointSet ): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        for i in range(1, pointSet.shape[2]):
            pointSet[:,i] = self.observe( pointSet[:,i] )
        return pointSet

    def saveStep( self ):
        self.__stateEstimateList.append( self.getStateEstimate() )
        self.__obsEstimateList.append( self.getObsEstimate() )

    def step( self, measurement ):
        
        self.setSigmaPoints()
        self.__sigmaPoints = self.transformSet(self.__sigmaPoints)
        self.__obsForecast = self.observeSet(self.__sigmaPoints)
        
        self.__stateForecast = self.__sigmaPoints @ transpose(self.__weights)
        self.__stateCov = (self.__weights*(self.__sigmaPoints-self.__stateForecast)) @ transpose(self.__sigmaPoints-self.__stateForecast)
        self.__obsEstimate = self.__obsForecast @ transpose(self.__weights)
        self.__obsCov = (self.__weights*(self.__obsForecast-self.__obsEstimate)) @ transpose(self.__obsForecast-self.__obsEstimate) + self.__sigmaHk
        self.__crossCov = (self.__weights*(self.__sigmaPoints-self.__selfEstimate)) @ transpose(self.__obsForecast-self.__obsEstimate)

        self.__gain = self.__crossCov * inverse(self.__obsCov)

        self.__selfEstimate += self.__gain @ ( measurement - self.__obsEstimate )

        self.saveStep()
        
    def filter( self, signal, initialGuess = "default", verbose = False ): # signal is a list of numpy arrays
        
        if initialGuess == "default":
            self.__setStateEstimate( np.zeros((self.stateDim,1)) )
        else:
            self.__setStateEstimate( initialGuess )
            
        N = len(signal)
        for k in range(N):
            self.step( signal[:,k] )
            
            if verbose:
                print("    sigma-points : ",self.getSigmaPoints() )
                print("    uncertainty matrix : ",self.getSigmaHk() )
                print("    Kalman gain : ",self.getGain() )
                print("    last state estimate : ",self.getStateEstimate() )
                print("    associated predicted measure : ", self.getObsEstimate(), " ; real measure : ", signal[:,i] )
                print("    predicted/obsereved measure gap : ",np.abs( self.getObsEstimate() - signal[:,i] ) )
