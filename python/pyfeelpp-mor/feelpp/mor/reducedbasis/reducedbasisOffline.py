from .reducedbasis import *

class reducedbasisOffline(reducedbasis):
    """ Reduced basis for stationnary problem
    """


    def __init__(self, Aq, Fq, model, mubar, output_names=None, use_dual_norm=False):
        """Initializes the class

        Args:
            Aq (list): affine decomposition of the bilinear form (Aq = [A1, ..., AqA])
            Fq (list): affine decomposition of the rhs AND the outputs Fq[0] = [F1, ..., FqF], F[i] = [Si1, ..., SiQsi]
            model (feelmm.mor._toolboxmor): model of toolboxmor, initialized
            mubar (parameterSpaceElement): parameter mubar
            outputs_names (list): list of the names of the diffetent outputs, in the order given by the environment
        """

        super().__init__(model)

        self.setMubar(mubar)        # default parameter

        self.Aq = Aq                # affine decomposition of A
        self.Fq = Fq[0]             # affine decomposition of F
        self.Lkq= Fq[1:]            # affine decompositions of Lk

        self.Qa = len(Aq)                           # size of the decomposition of A
        self.Qf = len(Fq[0])                        # size of the decomposition of F
        self.QLk= [len(Lkq) for Lkq in self.Lkq ]   # sizes of the decompositions of Lk
        self.N_output = len(self.QLk)               # number of outputs
        if output_names is None:
            self.output_names = [f"output{i}" for i in range(self.N_output)]
        else:
            self.output_names = output_names        # names of the outputs


        self.SS = np.zeros((self.Qf, self.Qf))      # SS[p,p_] = (S^p, S^p_)

        self.NN = Aq[0].size[0]     # size of the FE problem

        A_tmp = self.assembleA(self.betaA_bar[0])
        AT_tmp = A_tmp.copy()
        AT_tmp.transpose()

        self.scal = 0.5*(A_tmp + AT_tmp)        # scalar product for the energy norm
        self.Abar = self.assembleA(self.betaA_bar[0])


        # KSP to solve
        self.KSP_TYPE = PETSc.KSP.Type.GMRES
        self.PC_TYPE = PETSc.PC.Type.GAMG

        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_SELF)
        self.ksp.setType(self.KSP_TYPE)
        self.reshist = {}

        def monitor(ksp, its, rnorm):
            self.reshist[its] = rnorm
            if self.worldComm.isMasterRank():
                print("[petsc4py] Iteration {} Residual norm: {}".format(its,rnorm))

        def monitor_nverbose(ksp, its, rnorm):
            self.reshist[its] = rnorm

        self.monitorFunc = [monitor_nverbose, monitor]

        self.ksp.setMonitor(monitor_nverbose)
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)


        if not use_dual_norm:
            E = SLEPc.EPS()
            E.create()
            E.setOperators(self.Abar, self.scal)

            # for mubarbar:
            #   Compute the smallest eignevalue and keep it in memory for online stage TODO
            #   now it is only made with {mubarbar} = {mubar}

            E.setFromOptions()
            E.setWhichEigenpairs(E.Which.SMALLEST_MAGNITUDE)
            E.setDimensions(1)


            E.solve()

            nCv = E.getConverged()
            if self.worldComm.isMasterRank():
                print(f"[slepc4py] number of (smaller) eigenvalues computed : {nCv}")

            self.alphaMubar = E.getEigenvalue(0).real
            if self.worldComm.isMasterRank():
                print(f"[reducedbasis] Constant of continuity : {self.alphaMubar}")
        

            betaA_bar_np = np.array(self.betaA_bar[0])

            def alphaLB(mu):
                # From a parameter
                betaMu = self.model.computeBetaQm(mu)[0][0]
                return self.alphaMubar * np.min( betaMu / betaA_bar_np )

            def alphaLB_(betaA):
                # From a decomposition
                return self.alphaMubar * np.min( betaA / betaA_bar_np )


            E.setFromOptions()
            E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
            E.setDimensions(1)
            E.solve()

            nCv = E.getConverged()
            if self.worldComm.isMasterRank():
                print(f"[slepc4py] number of eigenvalues computed : {nCv}")

            self.gammaMubar = E.getEigenvalue(0).real
            if self.worldComm.isMasterRank():
                print(f"[reducedbasis] Constant of coercivity : {self.gammaMubar}")


            def gammaUB(mu):
                # From a parameter
                betaMu = self.model.computeBetaQm(mu)[0][0]
                return self.gammaMubar * np.max( betaMu / betaA_bar_np )

            def gammaUB_(betaA):
                # From a decomposition
                return self.gammaMubar * np.max( betaA / betaA_bar_np )

        else:
            alphaLB = lambda mu: 1
            alphaLB_ = lambda beta: 1
            gammaUB = lambda mu: 1
            gammaUB_ = lambda beta: 1

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_
        self.gammaUB = gammaUB
        self.gammaUB_ = gammaUB_

        self.isInitialized = True

        if self.worldComm.isMasterRank():
            print("[reducedbasis] Offline rb initialized")


    def setVerbose(self, set=True):
        """Set verbose when system is solved

        Args:
            set (bool, optional): Defaults to True.
        """
        self.ksp.setMonitor(self.monitorFunc[set])


    def fromConstructed(self, W):
        """Create an object from a previously generated basis.
        The offline erros are not computed in this function

        Args:
            W (list of PETSc.Vec): reduced basis

        Returns:
            reducedBasis: an object with W as basis
        """
        rb = reducedbasis(self.Aq, self.Fq + self.Lkq, self.model, self.mubar, self.alphaLB)
        self.Z = W
        self.N = len(W)

        self.generateANq()
        self.generateLkNp()
        self.generateFNp()

        return rb


    """
    Handle Gram-Schmidt orthogonalization
    """
    def scalarX(self, u, v):
        """Return the ernegy scalar product associed to the prolem

        Args:
            u (PETSc.Vec): vector
            v (PETSC.Vec): second vector

        Returns:
            float: v.T @ A @ u
        """
        # return v.dot( u )   # v.T @ scal @ u
        return v.dot( self.scal * u )   # v.T @ scal @ u

    def normA(self, u):
        """Compute the energy norm of the given vector

        Args:
            u (PETSc.Vec): vector

        Returns:
            float: ||u||_X
        """
        return np.sqrt(self.scalarX(u, u))

    def orthonormalizeZ(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis
        (the optional argument is not needed)
        """
        self.Z[0] /= self.normA(self.Z[0])
        for n in range(1, len(self.Z)):
            s = self.Z[0].duplicate()
            s.set(0)
            for m in range(n):
                s += self.scalarX(self.Z[n], self.Z[m]) * self.Z[m]
            z_tmp = self.Z[n] - s
            self.Z[n] = z_tmp / self.normA(z_tmp)
        # if not (self.test_orth() == np.eye(self.N)).all() and nb < 2:
        if not (self.test_orth() ) and nb < 10:
            self.orthonormalizeZ(nb=nb+1)
        elif self.worldComm.isMasterRank():
            print(f"[reducedBasis] Gram-Schmidt orthonormalization done after {nb+1} step"+['','s'][nb>1])





    """
    Assembly of matrices and vectors
    """
    def assembleA(self, beta):
        """Assemble the matrix from a given vector of parameters
           A = sum_q beta[q]*Aq[q]

        Args:
            beta (list): list of parameters betaA

        Returns:
            PETSc.Mat: assemble matrix
        """
        assert( len(beta) == self.Qa ), f"Number of param ({len(beta)}) should be {self.Qa}"

        A = self.Aq[0].duplicate()
        for q in range(0, self.Qa):
            A += self.Aq[q] * beta[q]
        return A

    def assembleF(self, beta):
        """Assemble the rhs from a given vector of parameters
           F = sum_q beta[q]*Fq[q]

        Args:
            beta (list): list of parameters betaF

        Returns:
            PETSc.Vec: assemble vector
        """
        assert( len(beta) == self.Qf ), f"Number of param ({len(beta)}) should be {self.Qf}"

        F = self.Fq[0].duplicate()
        for q in range(0, self.Qf):
            F += self.Fq[q] * beta[q]
        return F

    def assembleLk(self, k, beta):
        """Assemble the rhs from a given vector of parameters
           Lk = sum_q beta[q]*Lkq[q]

        Args:
            k (int): number of the i-th output
            beta (list): list of parameters betaF

        Returns:
            PETSc.Vec: assemble vector
        """

        assert( len(beta) == self.QLk[k]), f"Number of param ({len(beta)}) should be {self.QLk[k]}"

        L = self.Lkq[k][0].duplicate()
        for p in range(self.QLk[k]):
            L += self.Lkq[k][p] * beta[p]
        return L



    """
    Offline generation of reduced basis
    """

    def solveFE(self, mat, rhs):
        """Solve the finite elements problem max * X = rhs

        Args:
            mat (PETSc.Mat): matrix to inverse
            rhs (PETSc.Vec): right hand side of the equation

        Returns:
            PETSc.Vec: assemble vector
        """
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)

        self.reshist = {}

        self.ksp.setOperators(mat)
        self.ksp.setConvergenceHistory()
        sol = rhs.duplicate()
        self.ksp.solve(rhs, sol)

        return sol


    def generateZ(self, mus, orth=True):
        """Generates the base matrix Z (of shape (NN,N))

        Args:
            mus (list of ParameterSpaceElement): list of N parameters to evalutate offline
            orth (bool, optional): orthonormalize or not the reduced basis. Defaults to True.
        """
        self.N = len(mus)
        self.Z = []
        v = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        v.setSizes(self.NN)
        v.setUp()
        v.set(0)
        v.assemble()


        for i,mu in enumerate(tqdm(mus,desc=f"[reducedBasis] Offline generation of the basis", ascii=False, ncols=120)):
            beta = self.model.computeBetaQm(mu)
            A = self.assembleA(beta[0][0])
            F = self.assembleF(beta[1][0][0])

            sol = self.solveFE(A, F)
            self.Z.append(sol)
            # print(i, self.Z[-1].min(), self.Z[-1].max())
            # print(self.reshist)


        if orth:
            self.orthonormalizeZ()


    def test_orth(self):
        """Tests is the matrix Z is orthonormal
            Computes the matrix of scalar products of vectors of the reduced basis
            The returned matrix should be equal to the identity

        Returns:
            np.ndarray: ((xi_i,xi_j)_X)_{0<=i<N,0<=j<N}
        """
        sc = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                sc[i,j] = np.round(self.scalarX(self.Z[i], self.Z[j]), decimals=3)
                if np.round(np.round(self.scalarX(self.Z[i], self.Z[j]), decimals=3)) != [0,1][i==j]:
                    return False
        # print("sc",sc)
        # print("TESTORTH",np.linalg.norm(sc - np.eye(self.N)))
        # for i in range(self.N):
        #     for j in range(self.N):
        #         if np.round(np.round(self.scalarX(self.Z[i], self.Z[j]), decimals=3)) != [0,1][i==j]:
        #             return False
        return True


    def generateANq(self):
        """Generate the reduced matrices ANq
        """
        self.ANq = np.zeros((self.Qa, self.N, self.N))
        for i, ksi in enumerate(self.Z):
            for j, ksi_ in enumerate(self.Z):
                for q in range(self.Qa):
                    self.ANq[q,i,j] = ksi_.dot( self.Aq[q] * ksi )

    def expandANq(self):
        """Expand the reduced matrices ANq
        """
        self.ANq = np.concatenate( (self.ANq, np.zeros((self.Qa, self.N-1, 1))), axis=2)  # N has already been increased
        self.ANq = np.concatenate( (self.ANq, np.zeros((self.Qa, 1, self.N))), axis=1)
        for q in range(self.Qa):
            for i, ksi in enumerate(self.Z):
                self.ANq[q, -1, i] = ksi.dot( self.Aq[q] * self.Z[-1] )
                self.ANq[q, i, -1] = self.Z[-1].dot( self.Aq[q] * ksi )


    def generateFNp(self):
        """Generate the reduced vectors FNp
        """
        self.FNp = np.zeros((self.Qf, self.N))
        for i, ksi in enumerate(self.Z):
            for q in range(self.Qf):
                self.FNp[q,i] = self.Fq[q].dot(ksi)

    def expandFNp(self):
        self.FNp = np.concatenate( (self.FNp, np.zeros((self.Qf,1))), axis=1)
        for q in range(self.Qf):
            self.FNp[q, -1] = self.Fq[q].dot(self.Z[-1])

    def generateLkNp(self):
        """Generate the reduced vectors LkNp
        """
        self.LkNp = [[] for k in range(self.N_output)]
        for k in range(self.N_output):
            self.LkNp[k] = np.zeros((self.QLk[k], self.N))
            for i, ksi in enumerate(self.Z):
                for p in range(self.QLk[k]):
                    self.LkNp[k][p,i] = self.Lkq[k][p].dot(ksi)

    def expandLkNp(self):
        """Expand the reduced vectors LkNp
        """
        for k in range(self.N_output):
            self.LkNp[k] = np.concatenate( (self.LkNp[k], np.zeros((self.QLk[k],1))), axis=1)
            for q in range(self.QLk[k]):
                self.LkNp[k][q, -1] = self.Lkq[k][q].dot(self.Z[-1])



    def computeOfflineReducedBasis(self, mus, orth=True):
        """Computes the reduced basis and reduces matrices from a set of parameters

        Args:
            mus (list of ParameterSpaceElement) : list of parameters
            orth (bool, optional): orthonormalize or not the reduced basis. Defaults to True.
        """
        self.generateZ(mus,orth)
        self.generateANq()
        self.generateLkNp()
        self.generateFNp()

    """
    Finite elements resolution
    """
    def getSolutionsFE(self, mu, beta=None, k=-1):
        """Computes the finite element solution (big problem)

        Args:
            mu (ParameterSpaceElement) : parameter
            beta (list, optional) : coefficients of the decomposition, if they have already been computed
            k (int, optional) : index of the output to be computed, if -1 the compliant output is computed

        Returns:
            tuple (PETSc.Vec,float) : (u, s)
        """
        if beta is None:
            beta = self.model.computeBetaQm(mu)
        A_mu = self.assembleA(beta[0][0])
        F_mu = self.assembleF(beta[1][0][0])

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.reshist = {}

        self.ksp.setOperators(A_mu)
        self.ksp.setConvergenceHistory()
        sol = F_mu.duplicate()
        self.ksp.solve(F_mu, sol)

        if k == -1:
            l = F_mu
        else:
            if 0 <= k and k < self.N_output:
                l = self.assembleLk(k, beta[1][k+1][0])

        return sol, sol.dot(l)


    """
    Error computation
    """
    def computeOfflineErrorRhs(self):
        """Compute the offline errors associated to right hand side, independant of N
        """
        self.Sp = []

        # compute solutions of scal * Sp = Fp
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)

        for Fp in self.Fq:
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(Fp, sol)
            self.Sp.append(sol)

        # compute matrix of scalar products
        for p,Fp in enumerate(self.Sp):
            for p_,Fp_ in enumerate(self.Sp):
                self.SS[p,p_] = self.scalarX(Fp, Fp_)

    def computeOfflineError(self):
        """Compute offline errors associated to the reduced basis, dependant of N
        """
        self.Lnq = {}

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)
        for n, ksi in enumerate(self.Z):
            for q, Aq in enumerate(self.Aq):
                sol = self.Fq[0].duplicate()
                self.reshist = {}
                self.ksp.solve(-Aq*ksi, sol)
                self.Lnq[n,q] = sol.copy()

        self.SL = np.zeros((self.Qa, self.Qf, self.N))
        self.LL = np.zeros((self.Qa, self.N, self.Qa, self.N))
        for q in range(self.Qa):
            for n in range(self.N):
                for p in range(self.Qf):
                    self.SL[q,p,n] = self.scalarX(self.Sp[p], self.Lnq[n,q])

                for n_ in range(self.N):
                    for q_ in range(self.Qa):
                        self.LL[q,n,q_,n_] = self.scalarX(self.Lnq[n,q], self.Lnq[n_,q_])

    def expandOffline(self):
        """Add errors to values computed in previous steps.
           Before this function is called, the last column of Z must be computed
        """
        # self.Sp and self.SS are independant of N, so they don't change
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)

        self.SL = np.concatenate( ( self.SL, np.zeros( (self.Qa, self.Qf, 1) ) ), axis=2 )
        self.LL = np.concatenate( ( self.LL, np.zeros( (self.Qa, 1, self.Qa, self.N) ) ), axis=1 )
        self.LL = np.concatenate( ( self.LL, np.zeros( (self.Qa, self.N+1, self.Qa, 1) ) ), axis=3 )

        for q, Aq in enumerate(self.Aq):
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(-Aq*self.Z[-1], sol)

            self.Lnq[self.N, q] = sol

        for q in range(self.Qa):
            for p in range(self.Qf):
                self.SL[q,p,-1] = self.scalarX(self.Sp[p], self.Lnq[self.N,q])

            for n in range(self.N+1):
                for q_ in range(self.Qa):
                    self.LL[q,n,q_,-1] = self.scalarX( self.Lnq[n,q], self.Lnq[self.N,q_] )
                    self.LL[q,-1,q_,n] = self.scalarX( self.Lnq[self.N,q], self.Lnq[n,q_] )



    def computeDirectError(self, mu, precalc=None):
        """compute a posteriori error using a direct method (costly), from a parameter mu

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already\
                been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function

        Returns:
            float: a posteriori error
        """

        if precalc is None:
            beta_ = self.model.computeBetaQm(mu)
            betaA = beta_[0][0]
            betaF = beta_[1][0][0]
            bataLk= beta_[1][1:]
            AN_mu = self.assembleAN(betaA)
            FN_mu = self.assembleFN(betaF)

            uN = np.linalg.solve(AN_mu, FN_mu)
        else:
            betaA = precalc["betaA"]
            betaF = precalc["betaF"]
            betaLk= precalc["betaLk"]
            uN = precalc["uN"]

        A_mu = self.assembleA(betaA)
        F_mu = self.assembleF(betaF)

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.scal)
        self.reshist = {}

        E = self.Fq[0].duplicate()
        u_proj = self.projFE(uN)

        self.ksp.solve(F_mu - A_mu * u_proj, E)

        return self.normA(E)



    def compareSols(self, mu):
        """Compare solutions between reduced basis and finite element method

        Args:
            mu (ParameterSpaceElement): parameter

        Returns:
            float: norm || u_EF - u_RB ||_L2
        """
        uN,_ = self.getSolutions(mu)
        u, _ = self.getSolutionsFE(mu)

        u_proj = self.projFE(uN)

        return (u_proj - u).norm() / u.norm()

    def projFE(self, uN):
        """Project the reduced solution in FE space
        NB : this function is costly because it makes a loop on the big dimension NN
        NB : this function only works in sequential

        Args:
            uN (np.ndarray): reduced solution of size N

        Returns:
            PETSc.Vec: projection in FE space of size NN
        """
        u_proj = self.Fq[0].duplicate()
        u_proj.set(0)
        for i in range(self.NN):
            for j in range(len(uN)):
                u_proj[i] += self.Z[j][i] * uN[j]
        u_proj.assemble()
        return u_proj

    def RB_toNumpy(self):
        """Return the basis as a NumPy matrix (the quantities as copied)
        """
        RB = np.zeros((self.NN, self.N))
        for i, u in enumerate(self.Z):
            RB[:,i] = u[:]
        return RB





    """
    Greedy algorithm
    """
    def greedy(self, mu_0, Dmu, eps_tol=1e-6, Nmax=40):
        """Generate the basis, using the greedy algorithm

        Args:
            mu_0 (ParameterSpaceElement): first parameter for the basis
            Dmu (list of ParameterSpaceElement): test sample, to take arguments
            eps_tol (float, optional): tolerance. Defaults to 1e-6.
            Nmax (int, optional): Maximal size of the basis. Defaults to 40.

        Returns:
            list of ParameterSpaceElement: parameters for the basis
        """
        if self.worldComm.isMasterRank():
            print("[reducedBasis] Start greedy algorithm")
        self.N = 0
        S = []
        self.Z = []
        Delta = 1 + eps_tol

        self.DeltaMax = []

        self.computeOfflineErrorRhs()

        mu = mu_0

        betas = {}
        for i,mu in enumerate(tqdm(Dmu, desc=f"[reducedBasis] Computing betas", ascii=False, ncols=120)):
            betas[mu] = self.model.computeBetaQm(mu)

        # fig = go.Figure()

        # t_init = time.process_time()
        # ts = []

        while Delta > eps_tol and self.N < Nmax:

            S.append(mu)

            # self.computeOfflineReducedBasis(S)
            # self.computeOfflineError()

            if self.N == 0:
                self.computeOfflineReducedBasis(S)
                self.computeOfflineError()
            else:
                beta = betas[mu]
                A = self.assembleA(beta[0][0])
                F = self.assembleF(beta[1][0][0])

                sol = self.solveFE(A, F)
                self.Z.append(sol)
                self.orthonormalizeZ()

                self.expandOffline()

                self.N += 1

                # self.generateANq()
                # self.generateFNp()
                # self.generateLNp()
                self.expandANq()
                self.expandFNp()
                self.expandLkNp()


            mu_max = 0
            i_max = 0
            Delta_max = -float('inf')

            # t_loop = time.process_time()
            # t_sol = 0
            # t_i = 0
            # t_nrm = 0
            # t_eb = 0
            # t_tst = 0

            for i,mu_tmp in enumerate(tqdm(Dmu, desc=f"[reducedBasis] Greedy, step {self.N}", ascii=False, ncols=120)):
                beta = betas[mu_tmp]
                ANmu = self.assembleAN(beta[0][0])
                # t_sol0 = time.process_time()
                uN,_ = self.getSolutions(mu_tmp, beta=beta)
                # t_sol1 = time.process_time()
                norm_uMu = np.sqrt( uN.T @ ANmu @ uN )
                # t_nrm1 = time.process_time()

                # ti0 = time.process_time()
                precalc = {"betaA":beta[0][0], "betaF":beta[1][0][0], "uN":uN}
                # ti1 = time.process_time()

                # t_eb0 = time.process_time()
                D_en = self.computeEnergyBound(mu_tmp, precalc=precalc)
                Delta_tmp = D_en / norm_uMu
                # t_eb1 = time.process_time()

                # m.append(mu_tmp.parameterNamed('k_1'))
                # l.append(Delta_tmp)

                # t_tst0 = time.process_time()
                if Delta_tmp > Delta_max:
                    i_max = i
                    Delta_max = Delta_tmp
                    mu_max = mu_tmp
                # t_tst1 = time.process_time()

                # t_sol += t_sol1 - t_sol0
                # t_nrm += t_nrm1 - t_sol1
                # t_i += ti1 - ti0
                # t_eb += t_eb1 - t_eb0
                # t_tst += t_tst1 - t_tst0

            # t_end_loop = time.process_time()
            # ts.append((t_end_loop - t_loop))
            # print(ts[-1], ts[-1]/len(Dmu))
            # print(self.N, (t_sol)/len(Dmu), (t_nrm)/len(Dmu), (t_eb)/len(Dmu), (t_tst)/len(Dmu))
            Dmu.pop(i_max)
            Delta = Delta_max
            self.DeltaMax.append(Delta_max)
            mu = mu_max

            if self.worldComm.isMasterRank():
                print(f"[reducedBasis] Greedy algo, N={self.N}, Δ={Delta:e} (tol={eps_tol:e})".ljust(64), f"µ={mu}")
        if self.worldComm.isMasterRank() and self.N == Nmax:
            print("[reducedBasis] Greedy algo : warning, max size reached")

    #     fig.update_layout(
    #         # title="Energy for ",
    #         xaxis_title=r"$\mu$",
    #         yaxis_title=r"$\Delta_N$",
    #         legend_title=r"$N$"
    #     )
    #     fig.write_image("energy.png", scale=2)
    #     fig.write_html("energy.html")
        if self.worldComm.isMasterRank():
            print("[reducedBasis] End greedy algorithm")
        self.DeltaMax = np.array(self.DeltaMax)

        # t_end = time.process_time()
        # print(ts)
        # print(t_end - t_init)
        return S



    def generatePOD(self, Xi_train, output=-1, eps_tol=1e-6):
        """Generated the reduced basis using POD algorithm

        Args:
            Xi_train (list): set of train parameters
            output (int, optional): quantity of interest for the output. Defaults to -1.
            eps_tol (float, optional): tolerance for the RIC. Defaults to 1e-6.
        """

        if self.worldComm.isMasterRank():
            print("[reducedBasis] Start POD algorithm")
        self.N = 0
        self.Z = []
        Delta = 1 + eps_tol

        self.DeltaMax = []

        betas = {}
        solEF = {}
        for i,mu in enumerate(tqdm(Xi_train,desc=f"[reducedbasis] Computing offline solutions", ascii=False,ncols=120)):
            betas[mu] = self.model.computeBetaQm(mu)
            solEF[mu], _ = self.getSolutionsFE(mu, betas[mu])

        N_train = len(Xi_train)
        C = PETSc.Mat().create()
        C.setSizes([N_train, N_train])
        C.setFromOptions()
        C.setUp()

        for i,mu in enumerate(Xi_train):
            for j,nu in enumerate(Xi_train):
                C[i,j] = self.scalarX(solEF[mu], solEF[nu])
        C.assemble()

        E = SLEPc.EPS()
        E.create()
        E.setOperators(C)
        E.setFromOptions()
        E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
        E.setDimensions(N_train)

        E.solve()
        nCv = E.getConverged()
        if self.worldComm.isMasterRank():
            print(f"[slepsc4py] {nCv} eigenmodes computed (should be {N_train})")

        eigenval = np.zeros(nCv)
        for i in range(nCv):
            eigenval[i] = float(E.getEigenvalue(i).real)

        eigenvect = PETSc.Vec().create()
        eigenvect.setSizes(nCv)
        eigenvect.setFromOptions()
        eigenvect.setUp()

        while Delta > eps_tol:
            self.N += 1
            phi = solEF[Xi_train[0]].copy()
            phi.set(0)
            for i in range(nCv):
                E.getEigenvector(self.N-1, eigenvect)
                p = float(eigenvect[i]) * solEF[Xi_train[i]]
                phi += p
            self.Z.append(phi.copy())

            Delta = 1 - eigenval[:self.N].sum()/eigenval.sum()
            self.DeltaMax.append(Delta)

        self.orthonormalizeZ()
        self.generateANq()
        self.generateLkNp()
        self.generateFNp()
