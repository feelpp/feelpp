import numpy as np
from numpy import mean, sum, dot, linalg, zeros, ones, transpose
from numpy.linalg import inv, cholesky
from scipy.linalg import sqrtm

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
        self.P = np.eye(dim)

        self.SigPts = zeros([dim,2*dim+1])
        self.PreMeas = zeros([obs,2*dim+1])
        
        self.weights = ones(2*dim+1)*(1-w0)/(2*dim)
        self.weights[0] = w0
        self.covydefect = defect
        
        self.dynamics = dynamics
        self.observe = observe

    def readsignal(self, signal):
        self.signal = signal
        self.X = zeros([self.dim,max(signal.shape)-1]) # KEEPS TRACK OF THE STATES BEST ESTIMATE
        self.forecast = zeros([self.dim,max(signal.shape)]) # KEEPS TRACK OF THE FORECAST

    def inverse(M):          
        if M.size == 1:
            return 1/M
        else:
            return inv(M)
    
    def step(self):
        # COMPUTE SIGMA-POINTS
        self.P = sqrtm(self.P*self.dim/(1-self.weights[0]))
        self.SigPts[:,0] = self.Mx
        for i in range(1,self.dim+1):
            self.SigPts[:,i] = self.Mx + self.P[:,i-1]
        for i in range(self.dim+1,2*self.dim+1):
            self.SigPts[:,i] = self.Mx - self.P[:,i-self.dim-1]

        # TIME UPDATE
        for i in range(0,2*self.dim+1):
            self.SigPts[:,i] = self.dynamics(self.SigPts[:,i],self.Time,self.dt)
            self.PreMeas[:,i] = self.observe(self.SigPts[:,i])

        self.Time += 1

        self.Mx = self.SigPts @ transpose(self.weights)
        self.XF = self.Mx
        centerx = transpose([self.Mx]) @ [np.ones(2*self.dim+1)] * 1/(2*self.dim+1)
        self.Covx = (self.weights*(self.SigPts-centerx)) @ transpose(self.SigPts-centerx)
        self.My = self.PreMeas @ transpose(self.weights)
        centery = transpose([self.My]) @ [np.ones(2*self.dim+1)] * 1/(2*self.dim+1)
        self.Covy = (self.weights*(self.PreMeas-centery)) @ transpose(self.PreMeas-centery) + self.covydefect*np.eye(self.obs)
        self.XCov = (self.weights*(self.SigPts-centerx)) @ transpose(self.PreMeas-centery)

        # ANALYSIS STEP
        self.Kalman = self.XCov * self.inverse(self.Covy)
        self.Mx += self.Kalman @ (self.signal[self.Time]-self.My)
        self.P = self.Covx - self.Kalman @ np.transpose(self.XCov)
        
    def filter(self):
        for i in range(max(self.signal.shape)-1):
            self.X[:,i] = np.transpose(self.Mx)
            self.forecast[:,i] = np.transpose(self.XF)
            self.step(self)
            print(self.SigPts[0])
            print("Mx : ",self.Mx," My : ",self.My," measurement : ",self.signal[i])
            if np.abs(self.My-self.signal[i])/self.signal[i] < self.tol:
                return self.Mx
