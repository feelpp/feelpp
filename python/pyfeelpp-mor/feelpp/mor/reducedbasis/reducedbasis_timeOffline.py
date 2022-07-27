from .reducedbasis_time import *
from .reducedbasisOffline import *
import scipy.linalg as sl

class reducedbasisTimeOffline(reducedbasisOffline, reducedbasisTime):
    """
    Class for the offline part of the reduced basis method for the time-dependent problem
    """
    def __init__(self, *, Mr, **kwargs):
        """Initialize the reduced basis method for the time-dependent problem

        Args:
            Mr (list of PETScMat) : affine decomposition of the mass matrix M
            **kwargs              : keyword arguments for the reducedbasisTime class
        """
        warnings.warn("This class is still in construction. Correction from the previous version are in progress")
        super().__init__(**kwargs)

        self.Mr = Mr
        self.Qm = len(Mr)
        if self.Qm > 0:
            warnings.warn(f"The decomposition of the mass matrix should be of size 1, not {self.Qm}.\
                When assembleM will be called, only Mr[0] will be returned")

        self.Fkp : dict # size K*Qf : Fkp[k,p] <-> F^{k,p}
        self.Mnr : dict # size N*Qm : Mnr[k,p] <-> M^{n,r}


    def assembleM(self, beta):
        """Assemble the matrix from a given vector of parameters
           M = sum_r beta[r]*Mr[r]

        Args:
            beta (list): list of parameters betaM

        Returns:
            PETSc.Mat: assembled matrix
        """
        assert( len(beta) == self.Qm ), f"Number of param ({len(beta)}) should be {self.Qm}"

        # The decopmposition of the mass matrix is not correct yet
        # M = self.Mr[0].duplicate()
        # for r in range(0, self.Qm):
        #     M += self.Mr[r] * beta[r]
        return self.Mr[0]


    def computeOfflineReducedBasis(self, mus, orth=True):
        """Computes the reduced basis from a set of parameters

        Args:
            mus (list of ParameterSpaceElement): list of parameters
            orth (bool, optional): orthonormalize the basis. Defaults to True.
        """
        super().computeOfflineReducedBasis(mus, orth=orth)
        self.generateMNr()

    
    def generateMNr(self) -> None:
        """Generates the reduced matrices MNr
        """
        self.MNr = []
        for _ in range(self.Qm):
            self.MNr.append(np.zeros((self.N, self.N)))
        for i,u in enumerate(self.Z):
            for j,v in enumerate(self.Z):
                for r in range(self.Qm):
                    self.MNr[r][i,j] = v.dot( self.Mr[r] * u)


    def computeOfflineError(self, g):
        """Stores the offline data for error bound computation

        Args:
            g (np.ndarray): right-hand side time-dependant function
        """
        super().computeOfflineError()   # compute LL

        self.Fkp = {}
        self.Mnr = {}

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)
        sol = self.Fq[0].duplicate()

        for k in range(self.K):
            for p in range(self.Qf):
                rhs = self.Fq[p] * g((k+1)*self.dt)
                self.reshist = {}
                self.ksp.solve(rhs, sol)
                self.Fkp[k, p] = sol.copy()
        
        for n, ksi in enumerate(self.Z):
            for r, Mr in enumerate(self.Mr):
                self.reshist = {}
                self.ksp.solve( Mr * ksi, sol)
                self.Mnr[n, r] = sol.copy()
            
        self.FF = np.zeros((self.K, self.Qf, self.Qf))
        self.FM = np.zeros((self.K, self.Qf, self.Qm, self.N))
        self.FL = np.zeros((self.K, self.Qf, self.Qa, self.N))
        self.ML = np.zeros((self.Qm, self.N, self.Qa, self.N))
        self.MM = np.zeros((self.Qm, self.N, self.Qm, self.N))

        for k in range(self.K):
            for p in range(self.Qf):
                for p_ in range(self.Qf):
                    self.FF[k, p, p_] = self.scalarX(self.Fkp[k,p], self.Fkp[k,p_])
            
                for n in range(self.N):
                    for r in range(self.Qm):
                        self.FM[k,p,r,n] = self.scalarX(self.Fkp[k,p], self.Mnr[n,r])
                    for q in range(self.Qa):
                        self.FL[k,p,q,n] = self.scalarX(self.Lnq[n,q], self.Fkp[k,p])
        
        for n in range(self.N):
            for r in range(self.Qm):
                for n_ in range(self.N):
                    for q in range(self.Qa):
                        self.ML[r,n,q,n_] = self.scalarX(self.Lnq[n,q], self.Mnr[n,r])
                    for r_ in range(self.Qm):
                        self.MM[r,n,r_,n_] = self.scalarX(self.Mnr[n,r], self.Mnr[n_,r_])

    def expandOffline(self):
        super().expandOffline()
        # self.Fkp and self.FF are independant of N, so they don't change
        self.FM = np.concatenate( (self.FM,  np.zeros( (self.K, self.Qf, self.Qm, 1)   )), axis=3 )
        self.FL = np.concatenate( (self.FL,  np.zeros( (self.K, self.Qf, self.Qa, 1)   )), axis=3 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, 1, self.Qa, self.N)   )), axis=1 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, self.N+1, self.Qa, 1) )), axis=3 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, 1, self.Qm, self.N)   )), axis=1 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, self.N+1, self.Qm, 1) )), axis=3 )

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)

        for r, Mr in enumerate(self.Mr):
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(-Mr*self.Z[-1], sol)
            self.Mnr[self.N, r] = sol.copy()
        
        for k in range(self.K):
            for p in range(self.Qf):
                for r in range(self.Qm):
                    self.FM[k,p,r,-1] = self.scalarX(self.Fkp[k,p], self.Mnr[self.N,r])
                for q in range(self.Qa):
                    self.FL[k,p,q,-1] = self.scalarX(self.Lnq[self.N,q], self.Fkp[k,p])
        
        for r in range(self.Qm):
            for q in range(self.Qa):
                self.ML[r, -1, q, -1] = self.scalarX(self.Lnq[self.N,q], self.Mnr[self.N,r])
            for r_ in range(self.Qm):
                for n in range(self.N):
                    self.MM[r,n,r_,-1] = self.scalarX(self.Mnr[n,r], self.Mnr[self.N,r_])
                    self.MM[r,-1,r_,n] = self.scalarX(self.Mnr[self.N,r], self.Mnr[n,r_])


    """
    Offline generation of the basis
    """
    def generateBasis(self, musk, g=lambda k:1 if k==0 else 0, orth=True) -> None:
        """Generates the reduced basis matrix from different parameters and different instants

        Args:
            musk (dict): dict as {mu:[k0,k1,...], ...} where ki are sorted instants
        """

        warnings.warn("reducedbasis_time::reducedbasisTimeOffline::generateBasis has not yet been corrected or tested")

        self.N = taille_dict(musk)
        self.Z = []
        ind = 0

        for mu in musk:
            beta = self.model.computeBetaQm(mu)
            Amu = self.assembleA(beta[0][0])
            Fmu = self.assembleF(beta[1][0][0])
            Mmu = self.assembleM(beta[2][0])
            mat = Mmu + self.dt * Amu

            u = Fmu.duplicate()
            u.set(0)    # initial solution TODO

            mx = musk[mu][-1]
            cur = 0

            for k in range(mx):
                rhs = g(k) * self.dt * Fmu + Mmu * u
                self.ksp.setOperators(mat)
                self.ksp.setConvergenceHistory()
                sol = Fmu.duplicate()
                sol.set(0)
                self.ksp.solve(rhs, sol)
                u = sol.copy()

                if k+1 == musk[mu][cur]:
                    self.Z.append(sol.copy())
                    ind += 1
                    cur += 1

        if orth:
            self.orthonormalizeZ()
        self.generateANq()
        self.generateFNp()
        self.generateLkNp()
        self.generateMNr()


    """Functions dealing with the POD-greedy algorithm
    """

    def computeCorrelationMatrix(self, mu, g, to_numpy=True):
        """Computed the correlation matrix for POD-greedy algorithm

        Args:
            mu (ParameterSpaceElement): parameter used
            g (function): rhs function
            to_numpy (bool): if True, returns a numpy array, else a PETSc matrix

        Returns:
            np.ndarray (or PETSc.Mat): Correlation matrix C(mu)
        """
        eK = PETSc.Mat().create()
        eK.setSizes([self.NN, self.K])
        eK.setFromOptions()
        eK.setUp()

        uk = []

        betaA, betaF, betaM = self.model.computeBetaQm(mu)
            
        Amu = self.assembleA(betaA[0])
        Fmu = self.assembleF(betaF[0][0])
        Mmu = self.assembleM(betaM[0])

        mat = Mmu + self.dt * Amu
        self.ksp.setOperators(mat)

        u = self.Fq[0].duplicate()
        u.set(0)

        self.Z_to_matrix()
        ZT = self.Z_matrix.duplicate()
        ZT = ZT.transpose()
        ZZTX = self.Z_matrix * ZT * self.scal

        # build the matrix of error projections e^k(mu) = u^k(mu) - Z * ZT * Abar * u^k(mu)
        for k in range(self.K):

            rhs = g((k+1)*self.dt) * self.dt * Fmu + Mmu * u
            self.ksp.setConvergenceHistory()
            sol = self.Fq[0].duplicate()
            sol.set(0)
            self.ksp.solve(rhs, sol)
            
            eK[:,k] = sol - ZZTX * sol
            u = sol.copy()
            uk.append(sol.copy())
        
        eK.assemble()
        eKT = eK.copy()
        eKT = eKT.transpose()

        C = (eKT * self.Abar * eK) * 1./self.K
        if to_numpy:
            return C[:,:], uk
        else:
            return C, uk


    def computePODMode(self, mu, g, R=1):
        """compute the POD modes for a given parameter

        Args:
            mu (parameterSpaceElement): Paramter used
            g (function): function
            R (int, optional): Number of POD modes to compute. Defaults to 1.

        Returns:
            list: list of POD modes
        """
        C, uk = self.computeCorrelationMatrix(mu, g)
        values, vector = sl.eigh(C)
        ind = np.argsort(values)[::-1]
        res = []
        # TODO : add only 90% of the modes
        for r in range(R):
            psi_max = vector[ind[r]]
            POD = self.Z[0].duplicate()
            POD.set(0)
            for k in range(self.K):
                POD += float(psi_max[k]) * uk[k]
            res.append(POD)
        return res

    def greedyStep(self, xi_train, betas, g):
        mu_max = None
        i_max = 0
        Delta_max = -np.float('inf')
        for i,mu in enumerate(tqdm(xi_train,desc=f"[reducedBasis] POD-greedy, greey step", ascii=False, ncols=120)):

            uN = self.solveTime(mu, g)

            Delta = self.computeOnlineError(mu, betas, g, computeEnergyNorm=True)
            Delta_tmp = Delta / self.EnNorm[-1]

            if Delta_tmp > Delta_max:
                i_max = i
                mu_max = mu
                Delta_max = Delta_tmp
        return Delta_max, i_max, mu_max

    def generateBasisPODGreedy(self, mu0, mu_train, g, eps_tol=1e-6, R=1):
        """Run POD(t)-Greedy(µ) algorithm

        Args:
            mu0 (float): intitial parameter
            mu_train (list of ParameterSpaceElement): parameter train set
            g (function): function g
            eps_tol (float): critère d'arrêt.Default to 1e-6

        Returns:
            TO BE SPECIFIED
        """
        SN = []

        mu_star = mu0
        Delta_max = 1 + eps_tol

        betas = {}
        for i,mu in enumerate(tqdm(mu_train, desc=f"[reducedBasis] Computing betas", ascii=False, ncols=120)):
            betas[mu] = self.model.computeBetaQm(mu)

        self.computeOfflineReducedBasis([mu0])
        self.computeOfflineErrorRhs()
        self.computeOfflineError(g)

        while Delta_max > eps_tol:

            mu = mu_star
            SN.append(mu_star)

            # POD(t) step
            POD = self.computePODMode(mu, g, R=R)
            for ksi in POD:
                self.Z.append(ksi)
                self.expandOffline()
                self.N += 1

            self.generateANq()
            self.generateFNp()
            self.generateLkNp()
            self.generateMNr()


            # Greedy(µ) step
            Delta_max, i_star, mu_star = self.greedyStep(mu_train, betas, g)

            Delta = Delta_max
            mu_train.pop(i_star)
            print(f"[reducedbasis] POD-Greedy algorithm, N={self.N}, Δ={Delta} (tol={eps_tol})")

    def solveTimeForStudy(self, mu, g):
        """Computes both RB and FE solutions for a given parameter and a given time-dependent function

        Args:
            mu (parameterSpaceElement): parameter
            g (np.ndarray): function g

        Returns:
            tuple: results
        """

        warnings.warn("reducedbasis_time::reducedbasisTimeOffline::solveTimeForStudy has not yet been corrected or tested")
        beta = self.model.computeBetaQm(mu)
        ANmu = self.assembleAN(beta[0][0])
        FNmu = self.assembleFN(beta[1][0][0])
        MNmu = np.array(self.assembleMN(beta[2][0]))

        Amu = self.assembleA(beta[0][0])
        Fmu = self.assembleF(beta[1][0][0])
        Mmu = self.assembleM(beta[2][0])

        mat = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(mat)

        bigMat = Mmu + self.dt * Amu
        self.ksp.setOperators(bigMat)

        uN = np.zeros(self.N)
        u = self.Fq[0].duplicate()

        ones = self.Fq[0].duplicate()
        ones.set(1)
        onesN = np.ones(self.N)

        t = []
        sN = []
        s = []
        sDiff = []
        normN = []
        norm = []
        normDiff = []

        tmpN = np.zeros(self.K)
        tmp = np.zeros(self.K)
        tmpDiff = np.zeros(self.K)

        # résoudre dans le cas stationnaire au lieu de la boucle en temps TODO

        for k in range(1, self.K):
            t.append(k * self.dt)

            solN = sl.lu_solve(matLu, g(k*self.dt) * self.dt * FNmu + MNmu @ uN)
            uN = solN.copy()

            rhs = float(g[k]) * self.dt * Fmu + Mmu * u
            self.ksp.setConvergenceHistory()
            sol = self.Fq[0].duplicate()
            sol.set(0)
            self.ksp.solve(rhs, sol)
            u = sol.copy()

            sN.append(onesN @ MNmu @ uN)
            s.append(ones.dot(Mmu * u))
            sDiff.append(np.abs(s[-1] - sN[-1]))

            uplus = self.projFE(uN)
            diff = u - uplus
            
            tmpN[k] = uplus.dot( Amu * uplus )
            tmp[k] = u.dot( Amu * u )
            tmpDiff[k] = diff.dot(Amu * diff)

            normN.append( np.sqrt(uplus.dot(Mmu * uplus) + self.dt * tmpN.sum()) )
            norm.append( np.sqrt(u.dot(Mmu * u) + self.dt * tmp.sum()) )
            normDiff.append( np.sqrt(diff.dot(Mmu * diff) + self.dt * tmpDiff.sum()) )

        return t, sN, s, sDiff, normN, norm, normDiff
