from .reducedbasis_time import *
from .reducedbasisOffline import *

class reducedbasisTimeOffline(reducedbasisOffline, reducedbasisTime):
    """
    Class for the offline part of the reduced basis method for the time-dependent problem
    """
    def __init__(self, Aq, Fq, Mr, model, mubar, output_names=None, compute_lower_bound=True):
        warnings("This class is still in construction. Correction from the previous version are in progress")
        super().__init__(Aq, Fq, model, mubar, output_names, compute_lower_bound)

        self.Mr = Mr

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

        M = self.Mr[0].duplicate()
        for r in range(0, self.Qm):
            M += self.Aq[r] * beta[r]
        return M


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
                rhs = self.Fq[p] * g[k+1]
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
                    self.FF[k, p, p_] = self.scalarA(self.Fkp[k,p], self.Fkp[k,p_])
            
                for n in range(self.N):
                    for r in range(self.Qm):
                        self.FM[k,p,r,n] = self.scalarA(self.Fkp[k,p], self.Mnr[n,r])
                    for q in range(self.Qa):
                        self.FL[k,p,q,n] = self.scalarA(self.Lnq[n,q], self.Fkp[k,p])
        
        for n in range(self.N):
            for r in range(self.Qm):
                for n_ in range(self.N):
                    for q in range(self.Qa):
                        self.ML[r,n,q,n_] = self.scalarA(self.Lnq[n,q], self.Mnr[n,r])
                    for r_ in range(self.Qm):
                        self.MM[r,n,r_,n_] = self.scalarA(self.Mnr[n,r], self.Mnr[n_,r_])

    def expandOffline(self):
        super().expandOffline()
        # self.Fkp and self.FF are independant of N, so they don't change
        self.FM = np.concatenate( (self.FM, np.zeros(self.K, self.Qf, self.Qm, 1)), axis=3 )
        self.FL = np.concatenate( (self.FL, np.zeros(self.K, self.Qf, self.Qa, 1)), axis=3 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, 1, self.Qa, self.N) ) ), axis=1 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, self.N+1, self.Qa, 1) ) ), axis=3 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, 1, self.Qm, self.N) ) ), axis=1 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, self.N+1, self.Qm, 1) ) ), axis=3 )

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)

        for r, Mr in enumerate(self.Mr):
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(-Mr*self.Z[-1], sol)
            self.Mnr[self.N, r] = sol.copy()
        
        for k in range(self.K):
            for p in range(self.p):
                for r in range(self.Qm):
                    self.FM[k,p,r,-1] = self.scalarA(self.Fkp[k,p], self.Mnr[self.N,r])
                for q in range(self.Qa):
                    self.FL[k,p,q,-1] = self.scalarA(self.Lnq[self.N,q], self.Fkp[k,p])
        
        for r in range(self.Qm):
            for q in range(self.Qa):
                self.ML[r, -1, q, -1] = self.scalarA(self.Lnq[self.N,q], self.Mnr[self.N,r])
            for r_ in range(self.Qm):
                for n in range(self.N):
                    self.MM[-1,r,n,r_] = self.scalarA(self.Mnr[n,r], self.Mnr[self.N,r_])
                    self.MM[-1,r,n,r_] = self.scalarA(self.Mnr[self.N,r], self.Mnr[n,r_])


    """
    Offline generation of the basis
    """
    def generateBasis(self, musk, orth=True) -> None:
        """Generates the reduced basis matrix from different parameters and different instants

        Args:
            musk (dict): dict as {mu:[k0,k1,...], ...} where ki are sorted instants
        """

        warnings.warn("reducedbasis_time::reducedbasisTimeOffline::generateBasis has not yet been corrected or tested")

        self.N = taille_dict(musk)
        self.Z = []
        ind = 0

        def g(k): return 1 if k == 0 else 0
        # def g(k): return 1 - float(np.cos(k*self.dt))

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
        self.generateLNp()
        self.generateMNr()


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

            solN = matLu.solve(g[k] * self.dt * FNmu + MNmu @ uN)
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