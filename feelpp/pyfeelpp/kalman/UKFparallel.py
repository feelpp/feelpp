import numpy as np
from numpy import mean, sum, dot, linalg, zeros, ones, transpose
from numpy.linalg import inv, cholesky
from scipy.linalg import sqrtm
from mpi4py import MPI
from base import *

class Filter:

# CLASS INSTANCE SETTING IS MANDATORY : INITIALIZATION WITH STATE VECTOR DIMENSION

    def set(self, dim, obs, w0, dt, dynamics, observe, defect, initialguess, tol = 0.01):
        self.dim = dim
        self.obs = obs
        self.Time = 0
        self.dt = dt
        self.tol = tol
        
        self.Mx = initialguess               # mean value of x : state estimation
        self.Covx = zeros([dim,dim])
        self.My = ones(obs)                # mean value of y : signal estimation
        self.Covy = zeros([obs,obs])
        self.XCov = zeros([dim,obs])       # x y covariance
        self.XF = zeros(dim)               # current forecast

        self.Kalman = zeros([dim,obs])
        self.P = defect*np.eye(dim)

        self.SigPts = zeros([dim,2*dim+1], dtype=np.float64)
        self.PreMeas = zeros([obs,2*dim+1], dtype=np.float64)
        
        self.weights = ones(2*dim+1)*(1-w0)/(2*dim)
        self.weights[0] = w0
        self.covydefect = defect
        
        self.dynamics = dynamics
        self.observe = observe

        # MPI PARALLELIZATION INIT
        self.comm = MPI.COMM_WORLD
        self.nb_procs = self.comm.Get_size()
    
    def readsignal(self, signal):
        self.signal = signal
        self.X = zeros([self.dim,max(signal.shape)-1]) # KEEPS TRACK OF THE STATES BEST ESTIMATE
        self.forecast = zeros([self.dim,max(signal.shape)]) # KEEPS TRACK OF THE FORECAST

    def inverse(M):          
        if M.size == 1:
            return 1/M
        else:
            return inv(M)
    
    def step(self, mode):
        # COMPUTE SIGMA-POINTS
        self.P = sqrtm(self.P*self.dim/(1-self.weights[0]))
        self.SigPts[:,0] = self.Mx
        for i in range(1,self.dim+1):
            self.SigPts[:,i] = self.Mx + self.P[:,i-1]
        for i in range(self.dim+1,2*self.dim+1):
            self.SigPts[:,i] = self.Mx - self.P[:,i-self.dim-1]

        # TIME UPDATE

        # PARALLELIZED DYNAMICS
        rank = self.comm.Get_rank()
        
        counts = balancedpartition( 2*self.dim + 1, self.nb_procs )
        disp = displacements( counts )

        recvdata = np.zeros([2*self.dim + 1, counts[rank]], dtype=np.float64 )
        for ii in range( self.dim ): # SINCE SCATTER SENDS 1D ARRAY, PERFORMS LINEWISE SCATTER
            self.comm.Scatterv( [self.SigPts[ii], counts, disp, MPI.DOUBLE], recvdata[ii], root=0 )

        measdata = np.zeros([self.obs, counts[rank]])
        for i in range( 0, counts[rank] ): # USUAL POINTWISE PROPAGATION
            recvdata [:,i] = self.dynamics(recvdata[:,i],self.Time,self.dt)
            measdata[:,i] = self.observe(recvdata[:,i])

        sigptstmp = np.zeros([self.dim, counts[rank]], dtype=np.float64)
        premeastmp = np.zeros([self.obs, counts[rank]], dtype=np.float64)
        for ii in range( self.dim ):
            self.comm.Gatherv( [recvdata[ii], counts[rank]], [self.SigPts[ii], counts, disp, MPI.DOUBLE], root=0 )
        for ii in range( self.obs ):
            self.comm.Gatherv( [measdata[ii], counts[rank]], [self.PreMeas[ii], counts, disp, MPI.DOUBLE], root=0 )

        self.comm.Barrier()

        if (rank == 0):
            self.Time += 1

            self.Mx = self.SigPts @ transpose(self.weights)
            self.XF = self.Mx
            self.Covx = (self.weights*(self.SigPts-self.Mx)) @ transpose(self.SigPts-self.Mx)
            self.My = self.PreMeas @ transpose(self.weights)
            self.Covy = (self.weights*(self.PreMeas-self.My)) @ transpose(self.PreMeas-self.My) + self.covydefect*np.eye(self.obs)
            self.XCov = (self.weights*(self.SigPts-self.Mx)) @ transpose(self.PreMeas-self.My)

            # ANALYSIS STEP
            self.Kalman = self.XCov * self.inverse(self.Covy)
            
            if mode == "dynamic":
                self.Mx += self.Kalman @ (self.signal[self.Time]-self.My)
            elif mode == "static":
                self.Mx += self.Kalman @ (self.signal-self.My)
            
            self.P = self.Covx - self.Kalman @ np.transpose(self.XCov)
        
    def filter( self, measurement, maxiter = 1000, verbose = False, mode = "dynamic"):
        if mode == "dynamic":
            self.readsignal(self, measurement)
            for i in range(max(self.signal.shape)-1):
                self.X[:,i] = np.transpose(self.Mx)
                self.forecast[:,i] = np.transpose(self.XF)
                self.step(self, mode)
                if verbose:
                    print("    sigma-points : ",self.SigPts[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.Kalman)
                    print("    last state estimate : ",self.Mx)
                    print("    associated predicted measure : ",self.My," ; real measure : ",self.signal[i])
                    print("    relative measure error : ",np.abs(self.My-self.signal[i])/self.signal[i]," ; tolerance : ",self.tol)
                    if np.abs(self.My-self.signal[i])/self.signal[i] < self.tol or i < maxiter:
                        return self.Mx
                    
        elif mode == "static":
            i = 0   
            self.signal = measurement
            self.X = zeros([self.dim,maxiter]) # KEEPS TRACK OF THE STATES BEST ESTIMATE
            self.forecast = zeros([self.dim,maxiter]) # KEEPS TRACK OF THE FORECAST
            
            while np.abs(self.My-self.signal)/self.signal > self.tol and self.Time < maxiter:
                self.X[:,i] = np.transpose(self.Mx)
                self.forecast[:,i] = np.transpose(self.XF)
                self.step(self, mode)
                if verbose:
                    print("-------------------------------------------------------")
                    print("        TIME STEP",self.Time)
                    print(" ")
                    print("    sigma-points : ",self.SigPts[0])
                    print("    uncertainty matrix : ",self.P)
                    print("    Kalman gain : ",self.Kalman)
                    print("    last state estimate : ",self.Mx)
                    print("    associated predicted measure : ",self.My," ; real measure : ",self.signal)
                    print("    relative measure error : ",np.abs(self.My-self.signal)/self.signal," ; tolerance : ",self.tol)
                i += 1
                    
            return self.Mx
