import numpy
import scipy.linalg
from numpy.linalg import inv, cholesky
from base import *

class Filter:

    def __init__( self, dynamics = lambda x : x, observe = lambda x : x, defect = 1, stateDim = 1, obsDim = 1, w0 = 0.5 ):
        
        self.observe = observe
        self.transform = dynamics
        self.time = 0
        
        self.stateDim = stateDim
        self.obsDim = obsDim

        self._stateEstimateList = numpy.zeros((stateDim,0))
        self._obsEstimateList = numpy.zeros((obsDim,0))

        self._stateEstimate = numpy.ones((stateDim,1))    
        self._stateForecast = numpy.zeros((stateDim,1))
        self._obsEstimate = numpy.zeros((obsDim,1))

        self._obsForecast = numpy.zeros((obsDim,1+2*stateDim))
        self._stateCov = numpy.eye(stateDim)
        self._obsCov = numpy.eye(obsDim)
        self._crossCov = numpy.ones((stateDim,obsDim))
        self._sigmaHk = defect*numpy.eye(obsDim)
        self._sigmaAk = defect*numpy.eye(stateDim)
        self._gain = numpy.ones((stateDim,obsDim))
        
        self._sigmaScheme = defect*numpy.eye(stateDim)  # arbitrarily initialized with the observation defect parameter
        self._sigmaPoints = numpy.zeros((stateDim,2*stateDim+1))
        self._sigmaSigns = numpy.concatenate((numpy.zeros((stateDim,1)),numpy.eye(stateDim),-numpy.eye(stateDim)), axis=1)

        self._alpha, self._beta, self._kappa, self._lambda, self._factor = 0, 0, 0, 0, 0
        
        if 0: # method 1 - seems unstable
            self.setNumericalParameters( 1, 0, 3-stateDim )
            self._weightsMean = numpy.ones((1,2*stateDim+1))/(2*(stateDim+self._kappa))
            self._weightsMean[0] = self._kappa/(self.stateDim+self._kappa)
            self._weightsCov = self._weightsMean.copy()

        if 0: # method 2 - seems unstable
            self.setNumericalParameters( 1, 0, 3-stateDim )
            self._weightsMean = numpy.ones((1,2*stateDim+1))*self._lambda/(2*(stateDim+self._lambda))
            self._weightsCov = self._weightsMean.copy()
            self._weightsCov[0] += 1 - self._alpha**2 + self._beta
            self._factor = stateDim + self._lambda
            
        if 1: # method 3
            self._weightsMean = numpy.zeros((1,2*stateDim+1))
            self._weightsMean[0,0] = w0 # or any value < 1
            self._weightsMean[0,1:] = (1-self._weightsMean[0,0])/(2*stateDim)
            self._weightsCov = self._weightsMean
            self._factor = stateDim/(1-self._weightsMean[0,0])

        print("    Filter initialized with state dimension {} and observation dimension {}".format(stateDim, obsDim))

    def setNumericalParameters( self, a, b, k ):
        self._alpha, self._beta, self._kappa, self._lambda, self._factor = a, b, k, (self.stateDim + k)*a**2 - self.stateDim, self.stateDim + k

    def stepTime( self, dt ):
        self.time += dt

    def resetTime( self ):
        self.time = 0
        
    def setState( self, state ):
        self._stateEstimate = state
        
    def setSigmaHk( self, M ): # M is a numpy obsDim*obsDim matrix
        self._sigmaHk = M

    def setSigmaPoints( self ):
        self._sigmaScheme = numpy.sqrt(self._factor) * scipy.linalg.sqrtm( self._stateCov ).T @ self._sigmaSigns
        # numpy.linalg.cholesky( 3*self._stateCov ).T @ self._sigmaSigns
        self._sigmaPoints = self._stateEstimate + self._sigmaScheme

    def _setStateEstimate( self, value ):
        self._stateEstimate = value
        
    def getSigmaHk( self ):
        return self._sigmaHk
        
    def getSigmaPoints( self ):
        return self._sigmaPoints

    def getStateEstimate( self ):
        return self._stateEstimate.reshape(self.stateDim,1)

    def getStateForecast( self ):
        return self._stateForecast.reshape(self.stateDim,1)

    def getObsEstimate( self ): # this is the weighted mean of H(A(sigma-points))
        return self._obsEstimate.reshape(self.obsDim,1)

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
        return self.getStateEstimateList().shape[1]
        
    def transformSet( self, stateSet ): # stateSet is a numpy matrix representing points as columns ; to be parallelized
        self._sigmaPoints = numpy.apply_along_axis( self.transform, axis = 0, arr = stateSet, t = self.time )

    def observeSet( self, stateSet ): # stateSet is a numpy matrix representing points as columns ; to be parallelized
        self._obsForecast = numpy.apply_along_axis( self.observe, axis = 0, arr = stateSet, t = self.time )

    def _sigmaCov( self, x ):
        return numpy.outer( x, x )

    def _sigmaCrossCov( self, x, cut ):
        return numpy.outer( x[:cut], x[cut:] )

    def _weightedCov( self, x ):
        return numpy.inner( self._weightsCov, x )
    
    def saveStep( self ):
        self._stateEstimateList = numpy.hstack( ( self._stateEstimateList, self.getStateEstimate().copy() ) )
        self._obsEstimateList = numpy.hstack( (self._obsEstimateList,  self.getObsEstimate().copy() ) )

    def step( self, measurement ):
        
        self.setSigmaPoints()
        self.transformSet(self._sigmaPoints)
        self.observeSet(self._sigmaPoints)
        
        self._stateForecast = self._sigmaPoints @ self._weightsMean.T
        assert self._stateCov.dtype == "float64"

        self._stateCov = numpy.apply_along_axis(
            self._weightedCov,
            axis = 2,
            arr = numpy.apply_along_axis(
                self._sigmaCov,
                axis = 0,
                arr = self._sigmaPoints-self._stateForecast
            )
        ).reshape((self.stateDim,self.stateDim)) + self._sigmaAk

        assert self._stateCov.dtype == "float64"

        self._obsEstimate = self._obsForecast @ self._weightsMean.T
        
        self._obsCov =  numpy.apply_along_axis(
            self._weightedCov,
            axis = 2,
            arr = numpy.apply_along_axis(
                self._sigmaCov,
                axis = 0,
                arr = (self._obsForecast-self._obsEstimate).reshape(self.obsDim,1+2*self.stateDim)
            )
        ).reshape((self.obsDim,self.obsDim)) + self._sigmaHk

        assert self._obsCov.dtype == "float64"

        assert self._crossCov.dtype == "float64"
        
        self._crossCov = numpy.apply_along_axis(
            self._weightedCov,
            axis = 2,
            arr = numpy.apply_along_axis(
                self._sigmaCrossCov,
                axis = 0,
                arr = numpy.vstack(
                    (self._sigmaPoints - self._stateEstimate, self._obsForecast - self._obsEstimate)
                ).reshape(self.stateDim+self.obsDim, 1+2*self.stateDim),
                cut = self.stateDim
            )
        ).reshape((self.stateDim,self.obsDim))

        assert self._crossCov.dtype == "float64"

        self._gain = self._crossCov @ numpy.linalg.inv(self._obsCov)

        assert self._gain.dtype == "float64"

        self._stateCov -= self._gain @ self._obsCov @ self._gain.T
        
        self._stateEstimate += self._gain @ ( measurement.reshape((self.obsDim,1)) - self._obsEstimate )

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
