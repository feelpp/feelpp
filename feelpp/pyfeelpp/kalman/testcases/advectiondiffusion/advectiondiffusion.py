import numpy
from ukf import Filter
import matplotlib.pyplot

class Adv_Diff:

    def __init__( self, N = 101, dx = 0.01, dt = 0.005, mu = 0.01, v = 1, a = 1, b = 0 ):

        self.setDiscretization([N, dx, dt])
        self.setPhysicalParameters([mu, v])

        self._bc = [a, b]
        self._matrix = numpy.eye(N) - v*dt*(numpy.eye(N,N,1)+numpy.eye(N,N,-1))+mu*dt*(numpy.eye(N,N,1)+numpy.eye(N,N,-1)-2*numpy.eye(N))/(dx**2)
        self._state = numpy.zeros(N)
        self._state[[0,N-1]] = self._bc

    def step( self ):
        self._state = self._matrix @ self._state
        self.resetBC()

    def getDiscretization( self ):
        return [self._N, self._dx, self._dt]
    
    def getPhysicalParameters( self ):
        return [self._mu, self._v]

    def getBC( self ):
        return self._bc
    
    def getState( self ):
        return self._state

    def getMatrix( self ):
        return self._matrix

    def setDiscretization( self, triplet ):
        [self._N, self._dx, self._dt] = triplet

    def setPhysicalParameters( self, doublet ):
        [self._mu, self._v] = doublet

    def setBC( self, a, b ):
        self._bc = [a, b]
        
    def setState( self, x ):
        self._state = x
        
    def setMatrix( self, M ):
        self._matrix = M

    def resetBC( self ):
        self._state[[0,self._N-1]] = self._bc
        
    def plotState( self ):
        matplotlib.pyplot.plot( numpy.arange(0, 1+self._dx, self._dx), self.getState() )
        matplotlib.pyplot.show()


