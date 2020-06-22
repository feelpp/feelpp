import numpy
import scipy.linalg
from numpy.linalg import inv, cholesky
from base import *

class Filter:

    def __init__( self, dynamics = lambda x : x, observe = lambda x : x, defect = 1, stateDim = 1, obsDim = 1 ):
        self.observe = observe
        self.transform = dynamics
        
        self.stateDim = stateDim
        self.obsDim = obsDim

        self._stateEstimateList = []
        self._obsEstimateList = []
                
        self._stateEstimate = numpy.zeros((stateDim,1))
        self._stateForecast = numpy.zeros((stateDim,1))
        self._obsEstimate = numpy.zeros((obsDim,1))

        self._obsForecast = numpy.zeros((obsDim,1+2*stateDim))
        self._stateCov = numpy.eye(stateDim)
        self._obsCov = numpy.eye(obsDim)
        self._crossCov = numpy.ones((stateDim,obsDim))
        self._sigmaHk = defect*numpy.eye(obsDim)
        self._gain = numpy.ones((stateDim,obsDim))
        
        self._sigmaScheme = defect*numpy.eye(stateDim)  # arbitrarily initialized with the observation defect parameter
        self._sigmaPoints = numpy.zeros((stateDim,2*stateDim+1))
        self._sigmaSigns = numpy.concatenate((numpy.zeros((stateDim,1)),numpy.eye(stateDim),-numpy.eye(stateDim)), axis=1)
        
        self._weights = numpy.ones((1,2*stateDim+1))#*1/stateDim
#        self._weights[0] = 1-(stateDim-1)/stateDim # a negative weight causes issues with sqrt matrix
        self._weights3 = self._weights.copy()
        self._weights3.shape = (1,1,2*stateDim+1)
        
    def setSigmaHk( self, M ): # M is a numpy obsDim*obsDim matrix
        self._sigmaHk = M

    def setSigmaPoints( self ):
        self._sigmaScheme = scipy.linalg.sqrtm( 3*self._stateCov ).T @ self._sigmaSigns # numpy.transpose( numpy.linalg.cholesky( 3*self._stateCov ) ) @ self._sigmaSigns
#        print( self._stateCov )
        self._sigmaPoints = self._stateEstimate + self._sigmaScheme

    def _setStateEstimate( self, value ):
        self._stateEstimate = value
        
    def getSigmaHk( self ):
        return list(numpy.transpose(self._sigmaHk))
        
    def getSigmaPoints( self ):
        return self._sigmaPoints

    def getStateEstimate( self ):
        return self._stateEstimate

    def getStateForecast( self ):
        return self._stateForecast

    def getObsEstimate( self ): # this is the weighted mean of H(A(sigma-points))
        return self._obsEstimate

    def getObsForecast( self ): # this is the set of H(A(sigma-points))
        return self._obsForecast

    def getGain( self ):
        return self._gain

    def getWeights( self ):
        return self._weights

    def getStateEstimateList( self ):
        return self._stateEstimateList

    def getObsEstimateList( self ):
        return self._obsEstimateList

    def getNiter( self ):
        return len(self.getStateEstimateList())
        
    def transformSet( self, stateSet ): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        for i in range(stateSet.shape[1]):
            stateSet[:,i] = self.transform(stateSet[:,i])
        return stateSet

    def observeSet( self, stateSet ): # pointSet is a numpy matrix representing points as columns ; to be parallelized
        obsSet = np.zeros((self.obsDim,stateSet.shape[1]))
        for i in range(stateSet.shape[1]):
            obsSet[:,i] = self.observe( stateSet[:,i] )
        return obsSet

    def _sigmaCov( self, x ):
        return numpy.outer( x, x )
    
    def saveStep( self ):
        self._stateEstimateList.append( self.getStateEstimate().copy() )
        self._obsEstimateList.append( self.getObsEstimate().copy() )

    def step( self, measurement ):
        
        self.setSigmaPoints()
        self._sigmaPoints = self.transformSet(self._sigmaPoints)
        self._obsForecast = self.observeSet(self._sigmaPoints)
        
        self._stateForecast = self._sigmaPoints @ numpy.transpose(self._weights)
        self._stateCov = sum( ( self._weights3 * numpy.apply_along_axis( self._sigmaCov, axis = 0, arr = self._sigmaPoints - self._stateForecast ) ).T )

        self._obsEstimate = self._obsForecast @ numpy.transpose(self._weights)
        self._obsCov = sum( ( self._weights3 * numpy.apply_along_axis( self._sigmaCov, axis = 0, arr = self._obsForecast - self._obsEstimate ) ).T ) + self._sigmaHk
        # (self._weights*(self._obsForecast-self._obsEstimate)) @ numpy.transpose(self._obsForecast-self._obsEstimate) + self._sigmaHk
        self._crossCov = (self._weights*(self._sigmaPoints-self._stateEstimate)) @ numpy.transpose(self._obsForecast-self._obsEstimate)

        self._gain = self._crossCov @ numpy.linalg.inv(self._obsCov)
        print('===================================')
        print( self._stateEstimate.dtype, self._gain.dtype, measurement.dtype, self._obsEstimate.dtype )
        self._stateEstimate += self._gain @ ( numpy.matrix(measurement).T - self._obsEstimate )

        self.saveStep()
        
    def filter( self, signal, initialGuess = "default", verbose = False ): # signal is a list of numpy arrays
        
        if initialGuess == "default":
            self._setStateEstimate( numpy.zeros((self.stateDim,1)) )
        else:
            self._setStateEstimate( initialGuess )
            
        N = len(signal) # number of entries, whichever is their dimension
        for k in range(N):
            self.step( signal[k] )
            
            if verbose:
                print("================ step : ",k)
                print("    sigma-points : ",self.getSigmaPoints() )
                print("    uncertainty matrix : ",self.getSigmaHk() )
                print("    Kalman gain : ",self.getGain() )
                print("    last state estimate : ",self.getStateEstimate()[0] )
                print("    associated predicted measure : ", self.getObsEstimate()[0], " ; real measure : ", signal[k] )
                print("    predicted/obsereved measure gap : ",numpy.abs( self.getObsEstimate()[0] - signal[k] ) )
