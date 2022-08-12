from .reducedbasis_time import *
from .reducedbasisOffline import *
import scipy.linalg as sl
from tqdm import tqdm

def ric(val, N):
    total = val.sum()
    res = val[:N].sum()
    return res / total

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
        if self.Qm > 1:
            warnings.warn(f"The decomposition of the mass matrix should be of size 1, not {self.Qm}. When assembleM will be called, only Mr[0] will be returned")

        self.Mnr : dict # size N*Qm : Mnr[n,r] <-> M^{n,r}


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
        """Compute the reduced basis from a set of parameters

        Args:
            mus (list of ParameterSpaceElement): list of parameters
            orth (bool, optional): orthonormalize the basis. Defaults to True.
        """
        super().computeOfflineReducedBasis(mus, orth=orth)
        self.generateMNr()

    def generateMNr(self) -> None:
        """Generate the reduced matrices MNr
        """
        if self.Qm > 1:
            warnings.warn("Qm should have a decomposition of size 1")
        self.MNr = np.zeros((self.Qm, self.N, self.N))
        for i, ksi in enumerate(self.Z):
            for j, ksi_ in enumerate(self.Z):
                for r in range(self.Qm):
                    self.MNr[r,i,j] = ksi_.dot( self.Mr[r] * ksi)


    def computeOfflineError(self):
        """Store the offline data for error bound computation

        Args:
            g (function): right-hand side time-dependent function
        """
        super().computeOfflineError()   # compute LL

        self.Mnr = {}

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)
        sol = self.Fq[0].duplicate()

        for n, ksi in enumerate(self.Z):
            for r, Mr in enumerate(self.Mr):
                self.reshist = {}
                self.ksp.solve( -Mr * ksi, sol)
                self.Mnr[n, r] = sol.copy()

        self.FM = np.zeros((self.Qf, self.Qm, self.N))
        self.FL = np.zeros((self.Qf, self.Qa, self.N))
        self.ML = np.zeros((self.Qm, self.N, self.Qa, self.N))
        self.MM = np.zeros((self.Qm, self.N, self.Qm, self.N))

        for p in range(self.Qf):
            for n in range(self.N):
                for r in range(self.Qm):
                    self.FM[p,r,n] = self.scalarX(self.Sp[p], self.Mnr[n,r])
                for q in range(self.Qa):
                    self.FL[p,q,n] = self.scalarX(self.Lnq[n,q], self.Sp[p])

        for n in range(self.N):
            for r in range(self.Qm):
                for n_ in range(self.N):
                    for q in range(self.Qa):
                        self.ML[r,n,q,n_] = self.scalarX(self.Lnq[n,q], self.Mnr[n,r])
                    for r_ in range(self.Qm):
                        self.MM[r,n,r_,n_] = self.scalarX(self.Mnr[n,r], self.Mnr[n_,r_])

    def expandOffline(self):
        super().expandOffline()
        self.FM = np.concatenate( (self.FM,  np.zeros( (self.Qf, self.Qm, 1)   )), axis=2 )
        self.FL = np.concatenate( (self.FL,  np.zeros( (self.Qf, self.Qa, 1)   )), axis=2 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, 1, self.Qa, self.N)   )), axis=1 )
        self.ML = np.concatenate( ( self.ML, np.zeros( (self.Qm, self.N+1, self.Qa, 1) )), axis=3 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, 1, self.Qm, self.N)   )), axis=1 )
        self.MM = np.concatenate( ( self.MM, np.zeros( (self.Qm, self.N+1, self.Qm, 1) )), axis=3 )

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)

        for r, Mr in enumerate(self.Mr):
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(-Mr*self.Z[-1], sol)
            self.Mnr[self.N, r] = sol.copy()

        for p in range(self.Qf):
            for r in range(self.Qm):
                self.FM[p,r,-1] = self.scalarX(self.Sp[p], self.Mnr[self.N,r])
            for q in range(self.Qa):
                self.FL[p,q,-1] = self.scalarX(self.Lnq[self.N,q], self.Sp[p])

        for r in range(self.Qm):
            for q in range(self.Qa):
                self.ML[r, -1, q, -1] = self.scalarX(self.Lnq[self.N,q], self.Mnr[self.N,r])
            for r_ in range(self.Qm):
                for n in range(self.N+1):
                    self.MM[r,n,r_,-1] = self.scalarX(self.Mnr[n,r], self.Mnr[self.N,r_])
                    self.MM[r,-1,r_,n] = self.scalarX(self.Mnr[self.N,r], self.Mnr[n,r_])

    def computeDirectError(self, mu, beta, g, computeEnergyNorm=False):
        """Compute the error bound using the direct method

        Args:
            mu (float): regularization parameter
            beta (list): affine decomposition coefficients for parameter mu
            g (function): right-hand side time-dependent function
        """
        Amu = self.assembleA(beta[0][0])
        Fmu = self.assembleF(beta[1][0][0])
        Mmu = self.assembleM(beta[2][0])

        ANmu = self.assembleAN(beta[0][0])
        FNmu = self.assembleFN(beta[1][0][0])
        MNmu = self.assembleMN(beta[2][0])

        AZ = Amu * self.Z_matrix
        MZ = Mmu * self.Z_matrix
        matN = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(matN)
        alpm1 = 1 / np.sqrt(self.alphaLB(mu))

        uN = np.zeros(self.N)
        vN = PETSc.Vec().create()
        vN.setSizes(self.N)
        vN.setFromOptions()
        vN.setUp()
        vN_tmp = PETSc.Vec().create()
        vN_tmp.setSizes(self.N)
        vN_tmp.setFromOptions()
        vN_tmp.setUp()


        errDir = np.zeros(self.K)
        DeltaDir = np.zeros(self.K)
        aNorm = np.zeros(self.K)
        
        self.ksp.setOperators(self.scal)

        for k in range(1, self.K+1):
            gk = g(k * self.dt)
            uN_tmp = sl.lu_solve(matLu, gk * self.dt * FNmu + MNmu @ uN)

            self.ksp.setConvergenceHistory()
            sol = self.Fq[0].duplicate()
            vN.setValues(range(self.N), uN)
            vN_tmp.setValues(range(self.N), uN_tmp)

            t1 = float(gk) * Fmu
            v2 = t1.duplicate()
            MZ.mult( (vN_tmp-vN) / self.dt, v2) # v2 = MZ * (uN_tmp-uN) / self.dt
            v3 = t1.duplicate()
            AZ.mult(vN_tmp, v3)                 # v3 = AZ * uN_tmp

            self.ksp.solve(t1 - v2 - v3, sol)
            uN = uN_tmp.copy()

            aNorm[k-1] = uN.T @ ANmu @ uN
            errDir[k-1] = self.scalarX(sol, sol)
            DeltaDir[k-1] = alpm1 * np.sqrt(self.dt * errDir[:k].sum())

            if computeEnergyNorm:
                self.EnNorm[k-1] = np.sqrt(uN.T @ MNmu @ uN + self.dt * aNorm.sum())

        return DeltaDir[-1]

    """
    Offline generation of the basis
    """
    def generateBasis(self, musk, g=lambda k:1 if k==0 else 0, orth=True) -> None:
        """Generate the reduced basis matrix from different parameters and different instants

        Args:
            musk (dict): dict as {mu:[k0,k1,...], ...} where ki are sorted instants
            orth (bool, optional): orthonormalize the reduced basis. Defaults to True.
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
            u.set(0)    # initial condition TODO

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

    def computeCorrelationMatrix(self, mu, g, to_numpy=True, proj=True):
        """Compute the correlation matrix for POD-greedy algorithm

        Args:
            mu (ParameterSpaceElement): parameter used
            g (function): right-hand side time-dependent function
            to_numpy (bool): if True, returns a numpy array, else a PETSc matrix

        Returns:
            np.ndarray (or PETSc.Mat): Correlation matrix C(mu)
        """
        eK = PETSc.Mat().create()
        eK.setSizes((self.NN, self.K))
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
        u.set(0)     # initial condition TODO

        # build the matrix of error projections e^k(mu) = u^k(mu) - Z * ZT * scal * u^k(mu)
        if proj:
            self.Z_to_matrix()
            ZT = self.Z_matrix.copy()
            ZT = ZT.transpose()
            ZZTX = self.Z_matrix * ZT * self.scal

        for k in range(self.K):
            rhs = g((k+1)*self.dt) * self.dt * Fmu + Mmu * u
            self.ksp.setConvergenceHistory()
            sol = self.Fq[0].duplicate()
            sol.set(0)
            self.ksp.solve(rhs, sol)

            if proj:
                v = sol - ZZTX * sol
                # eK.setValuesBlocked(range(self.NN), k, sol - ZZTX * sol)
                for i in range(self.NN):
                    eK.setValue(i, k, v[i])
            else:
                # eK.setValuesBlocked(range(self.NN), k, sol)
                for i in range(self.NN):
                    eK.setValue(i, k, sol[i])

            u = sol.copy()
            uk.append(sol.copy())

        eK.assemble()
        eKT = eK.copy()
        eKT = eKT.transpose()

        C = (eKT * self.scal * eK) * 1./self.K      # C_{ij} = 1/K * ( u_i(µ), u_j(µ) )_X
        if to_numpy:
            return C[:,:], uk
        else:
            return C, uk


    def computePODMode(self, mu, g, R=1, delta=None, proj=True):
        """Compute the POD modes for a given parameter

        Args:
            mu (parameterSpaceElement): parameter used
            g (function): right-hand side time-dependent function
            R (int, optional): number of POD modes to compute. Defaults to 1.
            delta (float, optional): representativiness of the Nm first POD modes. If None, the R greatest POD modes are computed

        Returns:
            list: list of POD modes
        """
        C, uk = self.computeCorrelationMatrix(mu, g, proj=proj)
        values, vector = sl.eigh(C)
        ind = np.argsort(values)[::-1]
        res = []

        # Compute basis using R first POD modes
        if delta is None:
            for r in range(R):
                psi_max = vector[ind[r]]
                POD = self.Fq[0].duplicate()
                POD.set(0)
                for k in range(self.K):
                    POD += float(psi_max[k]) * uk[k]
                POD = POD * 1./np.sqrt(self.K)
                res.append(POD)
        # Compute basis using ric
        else:
            sum_eigen = values.sum()

            # Assert that the 3 largest eigenvalues have more than 90% of the energy
            if not ric(values[ind], 3) > 0.9:
                warnings.warn("The 3 largest eigenvalues have less than 90% of the energy")

            Nm = 0
            sum_delta = 0
            while sum_delta < delta * sum_eigen:
                psi_max = vector[ind[Nm]]
                POD = self.Fq[0].duplicate()
                POD.set(0)
                for k in range(self.K):
                    POD += float(psi_max[k]) * uk[k]
                POD = POD #* 1./np.sqrt(self.K)
                res.append(POD.copy())
                sum_delta += values[ind[Nm]]
                Nm += 1

        s = ["","s"][len(res)>1]
        print(f"[reducedbasis] POD-greedy, POD step : {len(res)} mode{s} computed")

        return res

    def greedyStep(self, xi_train, betas, g):
        """Run the greedy step to select the parameter

        Args:
            xi_train (list): list of parameters used for training
            betas (dict): dict containing the beta matrices for each parameter
            g (function): right-hand side time-dependent function

        Returns:
            tuple: Delta_max maximal error bound, i_max index of the parameter in xi_train, mu_max parameter selected
        """
        mu_max = None
        i_max = 0
        Delta_max = -np.float('inf')
        self.Z_to_matrix()
        for i,mu in enumerate(tqdm(xi_train,desc=f"[reducedBasis] POD-greedy, greedy step", ascii=False, ncols=120)):

            beta = betas[mu]

            # /!\ Z_to_matirx must be called before calling the next line
            # Delta = self.computeOnlineError(mu, beta, g, computeEnergyNorm=True)
            Delta = self.computeDirectError(mu, beta, g, computeEnergyNorm=True)
            Delta_tmp = Delta / self.EnNorm[-1]

            if Delta_tmp > Delta_max:
                i_max = i
                mu_max = mu
                Delta_max = Delta_tmp
        return Delta_max, i_max, mu_max

    def generateBasisPODGreedy(self, mu0, mu_train, g, eps_tol=1e-6, R=1, delta=None, doNotUseGreedy=False):
        """Run POD(t)-Greedy(µ) algorithm

        Args:
            mu0 (float): initial parameter
            mu_train (list of ParameterSpaceElement): parameter train set
            g (function): right-hand side time-dependent function
            eps_tol (float): stopping criterion.Default to 1e-6
            R (int, optional): number of POD modes to compute.
            delta (float, optional): representativiness of the Nm first POD modes

        Returns:
            list of parameters selected
        """
        SN = []
        self.Z = []

        mu_star = mu0
        Delta_max = 1 + eps_tol

        betas = {}
        for i,mu in enumerate(tqdm(mu_train, desc=f"[reducedBasis] Computing betas", ascii=False, ncols=120)):
            betas[mu] = self.model.computeBetaQm(mu)

        self.computeOfflineErrorRhs()

        while Delta_max > eps_tol:

            SN.append(mu_star)

            # POD(t) step
            if self.N == 0:
                POD = self.computePODMode(mu_star, g, R=R, delta=delta, proj=False)
                for ksi in POD:
                    self.Z.append(ksi)
                    self.N += 1
                self.computeOfflineError()
            else:
                POD = self.computePODMode(mu_star, g, R=R, delta=delta, proj=True)
                for ksi in POD:
                    self.Z.append(ksi)
                    self.expandOffline()
                    self.N += 1

            self.generateANq()
            self.generateFNp()
            self.generateLkNp()
            self.generateMNr()


            if doNotUseGreedy:
                i_star = np.random.randint(len(mu_train))
                mu_star = mu_train[i_star]
                Delta_max = Delta_max / 10
            else:
                # Greedy(µ) step
                Delta_max, i_star, mu_star = self.greedyStep(mu_train, betas, g)

            Delta = Delta_max
            mu_train.pop(i_star)
            print(f"[reducedbasis] POD-greedy algorithm, N={self.N}, Δ={Delta} (tol={eps_tol})")

        return SN

    def solveTimeForStudy(self, mu, g, k=-1):
        """Compute both RB and FE solutions for a given parameter and a given time-dependent function

        Args:
            mu (parameterSpaceElement): parameter used
            g (function): right-hand side time-dependent function
            k (int, optional): index of the output to be computed, if -1 the compliant output is computed. Default to -1

        Returns:
            tuple: results
        """

        beta = self.model.computeBetaQm(mu)

        ANmu = self.assembleAN(beta[0][0])
        FNmu = self.assembleFN(beta[1][0][0])
        MNmu = self.assembleMN(beta[2][0])

        Amu = self.assembleA(beta[0][0])
        Fmu = self.assembleF(beta[1][0][0])
        Mmu = self.assembleM(beta[2][0])

        mat = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(mat)

        bigMat = Mmu + self.dt * Amu
        self.ksp.setOperators(bigMat)

        uN = np.zeros(self.N)      # initial condition TODO
        u = self.Fq[0].duplicate() # initial condition TODO
        u.set(0)

        if k == -1:
            l = Fmu
            lN = FNmu
        else:
            if 0 <= k and k < self.N_output:
                l = self.assembleLk(k, beta[1][k+1][0])
                lN = self.assembleLkN(k, beta[1][k+1][0])
            else:
                raise ValueError(f"Output {k} not valid")

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
            gk = g(k * self.dt)

            solN = sl.lu_solve(matLu, gk * self.dt * FNmu + MNmu @ uN)
            uN = solN.copy()

            rhs = float(gk) * self.dt * Fmu + Mmu * u
            self.ksp.setConvergenceHistory()
            sol = self.Fq[0].duplicate()
            sol.set(0)
            self.ksp.solve(rhs, sol)
            u = sol.copy()

            sN.append( uN @ lN)
            s.append( u.dot(l) )
            sDiff.append(np.abs(s[-1] - sN[-1]))

            uplus = self.projFE(uN)
            diff = u - uplus

            tmpN[k] = uplus.dot( Amu * uplus )
            tmp[k] = u.dot( Amu * u )
            tmpDiff[k] = diff.dot( Amu * diff )

            normN.append( np.sqrt(uplus.dot(Mmu * uplus) + self.dt * tmpN.sum()) )
            norm.append( np.sqrt(u.dot(Mmu * u) + self.dt * tmp.sum()) )
            normDiff.append( np.sqrt(diff.dot(Mmu * diff) + self.dt * tmpDiff.sum()) )

        return t, sN, s, sDiff, normN, norm, normDiff

    """
    Checks functions, for debugging purposes
    """

    def checkDecomposition(self, mu, tol=1e-10):
        print("[reducedbasis] Checking decomposition with parameter µ={}".format(mu))
        beta = self.model.computeBetaQm(mu)

        ANmu = self.assembleAN(beta[0][0])
        FNmu = self.assembleFN(beta[1][0][0])
        MNmu = self.assembleMN(beta[2][0])

        Amu = self.assembleA(beta[0][0])
        Fmu = self.assembleF(beta[1][0][0])
        Mmu = self.assembleM(beta[2][0])

        vN = np.ones(self.N)
        vN_p = PETSc.Vec().createWithArray(vN)
        vN_p.setUp()

        v = self.Fq[0].duplicate()
        v.set(1)

        ZT = self.Z_matrix.copy()
        ZT = ZT.transpose()

        checksum = 0
        ncheck = 0

        # Check the matrix A
        # ANmu_p = Amu.ptap(self.Z_matrix)
        ANmu_p = ZT * Amu * self.Z_matrix
        check = vN.T @ ANmu @ vN
        check_p = vN.dot(ANmu_p * vN_p)
        print("    Check the matrix A", abs(check - check_p))
        checksum += abs(check - check_p); ncheck += 1

        # Check the matrix M
        # MNmu_p = Mmu.ptap(self.Z_matrix)
        MNmu_p = ZT * Mmu * self.Z_matrix
        check = vN.T @ MNmu @ vN
        check_p = vN.dot(MNmu_p * vN_p)
        print("    Check the matrix M", abs(check - check_p))
        checksum += abs(check - check_p); ncheck += 1

        # Check the matrix F
        # FNmu_p = Fmu.ptap(self.Z_matrix)
        FNmu_p = ZT * Fmu
        check = FNmu @ vN
        check_p = vN_p.dot( FNmu_p )
        print("    Check the vector F", abs(check - check_p))
        checksum += abs(check - check_p); ncheck += 1
        checksum += abs(check - check_p); ncheck += 1

        assert checksum / ncheck < tol, "Check failed"

    def checkError(self, mu, g, tol=1e-10):
        print("[reducedbasis] Checking error µ={}".format(mu))
        beta = self.model.computeBetaQm(mu)
        betaA = beta[0][0]
        betaF = beta[1][0][0]
        betaM = beta[2][0]

        Amu = self.assembleA(betaA)
        Fmu = self.assembleF(betaF)
        Mmu = self.assembleM(betaM)

        ANmu = self.assembleAN(betaA)
        FNmu = self.assembleFN(betaF)
        MNmu = self.assembleMN(betaM)

        AZ = Amu * self.Z_matrix
        MZ = Mmu * self.Z_matrix
        matN = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(matN)
        alpm1 = 1 / np.sqrt(self.alphaLB(mu))

        checksum = 0
        ncheck = 0

        uN_km1 = np.zeros(self.N)
        u_km1 = self.Fq[0].duplicate()
        u_km1.set(0)

        diff_p = PETSc.Vec().create()
        diff_p.setSizes(self.N)
        diff_p.setFromOptions()
        diff_p.setUp()
        u_p = diff_p.duplicate()
        rhs = self.Fq[0].duplicate()

        ONES = self.Fq[0].duplicate()
        ONES.set(1)

        pc = self.ksp.getPC()
        sol = self.Fq[0].duplicate()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)
        self.ksp.solve(Fmu, sol)

        sol_dec = self.Fq[0].duplicate()
        sol_dec.set(0)
        for p in range(self.Qf):
            sol_dec += betaF[p] * self.Sp[p]
        check = ONES.dot(sol)
        check_dec = ONES.dot(sol_dec)
        print("    Check the error on F", abs(check - check_dec))
        checksum += abs(check - check_dec); ncheck += 1

        for k in range(self.K):
            gk = g(k * self.dt)

            # computation of the solution
            uN_k = sl.lu_solve(matLu, gk * self.dt * FNmu + MNmu @ uN_km1)

            # computation of Mk
            diff_N = (uN_k - uN_km1) / self.dt
            diff_p.setValues(range(self.N), diff_N)
            (Mmu * self.Z_matrix).mult( diff_p, rhs )
            self.ksp.solve(-rhs, sol)
            Mk = sol.copy()

            # computation of Mk_dec
            Mk_dec = self.Fq[0].duplicate()
            Mk_dec.set(0)
            for r in range(self.Qm):
                for n in range(self.N):
                    Mk_dec += betaM[r] * self.Mnr[n,r] * diff_N[n]
            
            check = ONES.dot(Mk)
            check_dec = ONES.dot(Mk_dec)
            print(f"    Check the error on M{k}", abs(check - check_dec))
            checksum += abs(check - check_dec); ncheck += 1

            # computation of Lk
            u_p.setValues(range(self.N), uN_k)
            (Amu * self.Z_matrix).mult( u_p, rhs )
            self.ksp.solve(-rhs, sol)
            Lk = sol.copy()

            # computation of Lk_dec
            Lk_dec = self.Fq[0].duplicate()
            Lk_dec.set(0)
            for q in range(self.Qa):
                for n in range(self.N):
                    Lk_dec += betaA[q] * self.Lnq[n,q] * uN_k[n]

            check = ONES.dot(Lk)
            check_dec = ONES.dot(Lk_dec)
            print(f"    Check the error on L{k}", abs(check - check_dec))
            checksum += abs(check - check_dec); ncheck += 1

        assert checksum / ncheck < tol, "Check failed"
