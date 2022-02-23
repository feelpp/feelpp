"""reducedBasisPetsc :
    deals with reduced basis, using PETSc matrices
"""

import numpy as np
from petsc4py import PETSc
# import plotly.graph_objects as go
import sys, os
from tqdm import tqdm
from scipy.sparse.linalg import splu, spsolve



def convertToPetscMat(Aq):
    AqP = []
    for A in Aq:
        AqP.append(A.to_petsc().mat())
        if not AqP[-1].assembled:
            AqP[-1].assemble()
    return AqP

def convertToPetscVec(Fq):
    FqP = []
    for F in Fq:
        FqP.append(F.to_petsc().vec())
    return FqP


class reducedbasis():
    """ Reduced basis for stationnary problem
    """

    def __init__(self, Aq, Fq, model, mubar) -> None:
        """Initialise the object

        Args:
            Aq (list of PETSc.Mat): matrices Aq given by the affine decomposition
            Fq (list of PETSc.Vec): vectors Fq given by the decomposition of right-hand side
            model (ToolboxMor_{2|3}D): model DEIM used for the decomposition
            mubar (ParameterSpaceElement): parameter mu_bar for the enrgy norm
        """
        
        self.Aq = Aq    # of len Qa
        self.Fq = Fq    # of len Qf

        self.NN = Aq[0].size[0] if Aq is not None else 0
        self.Qa = len(Aq) if Aq is not None else 0
        self.Qf = len(Fq) if Fq is not None else 0
        self.N = 0

        self.model = model
        self.mubar = mubar

        self.Z      : list  # len N,  each vector size NN
        self.ANq    : list  # len Qa, each matrix size (N,N)
        self.FNp    : list  # len Qf, each vector size N

        self.Sp     : list  # len Qf, each vector of size NN
        self.Lnq    : dict  # size N*Qa : self.Lnq[n,q] <-> L^{n,q}

        self.SS = np.zeros((self.Qf, self.Qf))            # SS[p,p_] = (Sp,sp_)X
        self.SL     : np.ndarray    # shape (Qf, Qa, N)     SL[p,n,q] = (Sp,Lnq)X
        self.LL     : np.ndarray    # shape (Qa, N, Qa, N)  LL[n,p,n_,p_] = (Lnq,Ln_p_)X

        if self.mubar is not None:
            beta = self.model.computeBetaQm(self.mubar)
            self.betaA_bar = beta[0]
            self.betaF_bar = beta[1]
            A_tmp = self.assembleA(self.betaA_bar[0])
            AT_tmp = A_tmp.copy()
            AT_tmp.transpose()
            self.Abar = 0.5*(A_tmp + AT_tmp)

        # KSP to solve
        self.KSP_TYPE = PETSc.KSP.Type.GMRES
        self.PC_TYPE = PETSc.PC.Type.LU

        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_SELF)
        self.ksp.setType(self.KSP_TYPE)
        self.reshist = {}

        def monitor(ksp, its, rnorm):
            self.reshist[its] = rnorm
            print("[petsc4py] Iteration {} Residual norm: {}".format(its,rnorm))
        
        def monitor_nverbose(ksp, its, rnorm):
            self.reshist[its] = rnorm
        
        self.monitorFunc = [monitor_nverbose, monitor]

        self.ksp.setMonitor(monitor)
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)


        E = SLEPc.EPS()
        E.create()
        E.setOperators(self.Abar)

        E.setFromOptions()
        E.setWhichEigenpairs(E.Which.LARGEST_REAL)
        # E.setDimensions(5)

        E.solve()

        nCv = E.getConverged()
        if rank == 0:
            print(f"[slepc4py] number of eigenvalues computed : {nCv}")

        alphaMubar = E.getEigenvalue(0).real
        # alphaMubar = 1
        if rank == 0:
            print(f"[reducedbasis] Constant of continuity : {alphaMubar}")
        betaA_bar_np = np.array(self.betaA_bar)
        
        def alphaLB(mu):
            betaMu = self.model.computeBetaQm(self.mubar)[0][0]
            return alphaMubar * np.min( betaMu / betaA_bar_np )

        self.alphaLB = alphaLB

        ## For greedy memory
        self.DeltaMax = None

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
            W (list of PETSc.Vec): rediced basis

        Returns:
            reducedBasis: an object with W as basis
        """
        rb = reducedbasis(self.Aq, self.Fq, self.model, self.mubar, self.alphaLB)
        self.Z = W
        self.N = len(W)

        self.generateANq()
        self.generateLNp()
        self.generateFNp()

        return rb


    """
    Handle Gram-Schmidt orthogonalization
    """   
    def scalarA(self, u, v):
        """Return the ernegy scalar product associed to the matrix Abar

        Args:
            u (PETSc.Vec): vector
            v (PETSC.Vec): second vector

        Returns:
            float: v.T @ A @ u
        """
        # return v.dot( u )   # v.T @ Abar @ u
        return v.dot( self.Abar * u )   # v.T @ Abar @ u
    
    def normA(self, u):
        """Compute the energy norm of the given vector

        Args:
            u (PETSc.Vec): vector

        Returns:
            float: ||u||_X
        """
        return np.sqrt(self.scalarA(u, u))
    
    def orthonormalizeZ(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis
        (the optional argument is not needed)
        """
        self.Z[0] /= self.normA(self.Z[0])
        for n in range(1, len(self.Z)):
            s = self.Z[0].duplicate()
            s.set(0)
            for m in range(n):
                s += self.scalarA(self.Z[n], self.Z[m]) * self.Z[m]
            z_tmp = self.Z[n] - s
            self.Z[n] = z_tmp / self.normA(z_tmp)
        # if not (self.test_orth() == np.eye(self.N)).all() and nb < 2:
        if not (self.test_orth() ) and nb < 10:
            self.orthonormalizeZ(nb=nb+1)
        else:
            # pass
            print(f"[reducedBasis] Gram-Schmidt orthonormalization done after {nb+1} step"+['','s'][nb>0])


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

    def assembleAN(self, beta, size=None):
        """Assemble the reduced matrix from a given vector of parameters
           AN = sum_q beta[q]*ANq[q]

        Args:
            beta (list): list of parameters betaA
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            PETSc.Mat: assemble matrix
        """
        assert( len(beta) == self.Qa ), f"Number of param ({len(beta)}) should be {self.Qa}"

        AN = self.ANq[0].copy()
        AN.fill(0)
        for q in range(0, self.Qa):
            AN += self.ANq[q] * beta[q]
        if size is None:
            return AN
        else:
            return AN[:size,:size]

    def assembleFN(self, beta, size=None):
        """Assemble the rhs from a given vector of parameters
           FN = sum_q beta[q]*FNp[q]

        Args:
            beta (list): list of parameters betaF
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            PETSc.Vec: assemble vector
        """
        assert( len(beta) == self.Qf ), f"Number of param ({len(beta)}) should be {self.Qf}"

        FN = self.FNp[0].copy()
        FN.fill(0)
        for q in range(0, self.Qf):
            FN += self.FNp[q] * beta[q]
        if size is None:
            return FN
        else:
            return FN[:size]



    """
    Offline generation of reduced basis
    """

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
        i=0

        for mu in mus:
            beta = self.model.computeBetaQm(mu)
            A = self.assembleA(beta[0][0])
            F = self.assembleF(beta[1][0][0])
            
            pc = self.ksp.getPC()
            pc.setType(self.PC_TYPE)

            self.reshist = {}

            self.ksp.setOperators(A)
            self.ksp.setConvergenceHistory()
            sol = F.duplicate()
            self.ksp.solve(F, sol)
            self.Z.append(sol)
            # print(i, self.Z[-1].min(), self.Z[-1].max())
            # print(self.reshist)
            i += 1

        if orth:
            self.orthonormalizeZ()


    def test_orth(self):
        """Tests is the matrix Z is orthonormal
        #    Computes the matrix of scalar products of vectors of the reduced basis
        #    The returned matrix should be equal to the identity

        Returns:
            np.ndarray: ((xi_i,xi_j)_X)_{0<=i<N,0<=j<N}
        """
        sc = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                sc[i,j] = np.round(self.scalarA(self.Z[i], self.Z[j]), decimals=3)
                if np.round(np.round(self.scalarA(self.Z[i], self.Z[j]), decimals=3)) != [0,1][i==j]:
                    return False
        # print("sc",sc)
        # print("TESTORTH",np.linalg.norm(sc - np.eye(self.N)))
        # for i in range(self.N):
        #     for j in range(self.N):
        #         if np.round(np.round(self.scalarA(self.Z[i], self.Z[j]), decimals=3)) != [0,1][i==j]:
        #             return False
        return True


    def generateANq(self):
        """Generate the reduced matrices ANq
        """
        self.ANq = []
        for q in range(self.Qa):
            self.ANq.append(np.zeros((self.N,self.N)))
        for i,u in enumerate(self.Z):
            for j,v in enumerate(self.Z):
                for q in range(self.Qa):
                    self.ANq[q][i,j] = v.dot( self.Aq[q] * u )
                    # print(i, j, self.ANq[q][i,j]) # uT.A.u
    
    def generateLNp(self):
        """Generate the reduced vectors LNq
        """
        self.LNp = []
        for q in range(self.Qf):
            self.LNp.append(np.zeros((self.N)))
        for i,u in enumerate(self.Z):
            for q in range(self.Qf):
                self.LNp[q][i] = self.Fq[q].dot(u)
                # print(i, self.LNp[q][i])
    
    def generateFNp(self):
        """Generate the reduced vectors FNp
        """
        self.FNp = []
        for q in range(self.Qf):
           self.FNp.append(np.zeros((self.N)))
        for i,u in enumerate(self.Z):
            for q in range(self.Qf):
                self.FNp[q][i] = self.Fq[q].dot(u)


    def computeOfflineReducedBasis(self, mus, orth=True):
        """Computes the reduced basis and reduces matrices from a set of parameters

        Args:
            mus (list of ParameterSpaceElement) : list of parameters
            orth (bool, optional): orthonormalize or not the reduced basis. Defaults to True.
        """
        self.generateZ(mus,orth)
        self.generateANq()
        self.generateLNp()
        self.generateFNp()

    """
    Finite elements resolution
    """
    def getSolutionsFE(self, mu):
        """Computes the finite element solution (big problem)

        Args:
            mu (ParameterSpaceElement) : parameter

        Returns:
            tuple (PETSc.Vec,float) : (u, s)
        """
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

        return sol, sol.dot(F_mu)

    """
    Online computation
    """
    def getSolutions(self, mu, size=None):
        """Return solution uN and output sN

        Args:
            mu (ParameterSpaceElement) : parameter
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            tuple (np.ndarray,float) : (uN, sN)
        """
        beta = self.model.computeBetaQm(mu)
        A_mu = self.assembleAN(beta[0][0], size=size)
        F_mu = self.assembleFN(beta[1][0][0], size=size)

        sol = np.linalg.solve(A_mu, F_mu)

        return sol, sol @ F_mu

    """
    Error computation
    """
    def computeOfflineErrorRhs(self):
        """Compute the offline errors associated to right hand side, independant of N
        """
        self.Sp = []
        
        # compute solutions of Abar Sp = Fp
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)

        for Fp in self.Fq:
            sol = self.Fq[0].duplicate()
            self.reshist = {}
            self.ksp.solve(Fp, sol)
            self.Sp.append(sol)

        # compute matrix of scalar products
        for p,Fp in enumerate(self.Sp):
            for p_,Fp_ in enumerate(self.Sp):
                self.SS[p,p_] = self.scalarA(Fp, Fp_)

    def computeOfflineError(self):
        """Compute offline errors associated to the reduced basis, dependant of N
        """
        self.Lnq = {}

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)
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
                    self.SL[q,p,n] = self.scalarA(self.Sp[p], self.Lnq[n,q])

                for n_ in range(self.N):
                    for q_ in range(self.Qa):
                        self.LL[q,n,q_,n_] = self.scalarA(self.Lnq[n,q], self.Lnq[n_,q_])

    def expandOffline(self):
        """Add errors to values computed in previous steps.
           Before this function is called, the last column of Z must be computed
        """
        # self.Sp and self.SS are independant of N, so they don't change
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)

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
                self.SL[q,p,-1] = self.scalarA(self.Sp[p], self.Lnq[self.N,q])

            for n in range(self.N+1):
                for q_ in range(self.Qa):
                    self.LL[q,n,q_,-1] = self.scalarA( self.Lnq[n,q], self.Lnq[self.N,q_] )
                    self.LL[q,-1,q_,n] = self.scalarA( self.Lnq[self.N,q], self.Lnq[n,q_] )

        

            


    def computeOnlineError(self, mu, precalc=None):
        """compute a posteriori online error, from a parameter mu

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function

        Returns:
            float: a posteriori error SQUARED
        """

        if precalc is None:
            beta_ = self.model.computeBetaQm(mu)
            betaA = beta_[0][0]
            betaF = beta_[1][0][0]
            A_mu = self.assembleAN(betaA)
            F_mu = self.assembleFN(betaF)

            uN = np.linalg.solve(A_mu, F_mu)
        else:
            betaA = precalc["betaA"]
            betaF = precalc["betaF"]
            uN = precalc["uN"]

        s1 = betaF @ self.SS @ betaF    # beta_p*beta_p'*(Sp,Sp')
        s2 = 0
        s3 = 0

        # s1 = 0
        # for p in range(self.Qf):
        #     for p_ in range(self.Qf):
        #         s1 += betaA[p] * betaA[p_] * self.SS[p,p_]
        s2 = 0
        for p in range(self.Qf):
            for q in range(self.Qa):
                for n in range(self.N):
                    s2 += betaF[p] * betaA[q] * uN[n] * self.SL[q,p,n]

        s3 = 0
        for q in range(self.Qa):
            for q_ in range(self.Qa):
                for n in range(self.N):
                    for n_ in range(self.N):
                        s3 += betaA[q] * betaA[q_] * uN[n] * uN[n_] * self.LL[q,n,q_,n_]

        return s1 + 2*s2 + s3

    def computeDirectError(self, mu, precalc=None):
        """compute a posteriori error using a direct method (costly), from a parameter mu

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function

        Returns:
            float: a posteriori error
        """

        if precalc is None:
            beta_ = self.model.computeBetaQm(mu)
            betaA = beta_[0][0]
            betaF = beta_[1][0][0]
            AN_mu = self.assembleAN(betaA)
            FN_mu = self.assembleFN(betaF)

            uN = np.linalg.solve(AN_mu, FN_mu)
        else:
            betaA = precalc["betaA"]
            betaF = precalc["betaF"]
            uN = precalc["uN"]

        A_mu = self.assembleA(betaA)
        F_mu = self.assembleF(betaF)

        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)
        self.ksp.setOperators(self.Abar)
        self.reshist = {}
        
        E = self.Fq[0].duplicate()
        u_proj = self.projFE(uN)

        self.ksp.solve(F_mu - A_mu * u_proj, E)
        
        return self.normA(E)


    def computeEnergyBound(self, mu, precalc=None):
        """computes ernergy bound

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function

        Returns:
            float: energy error bound
        """
        normHatE = self.computeOnlineError(mu, precalc=precalc)
        alp = self.alphaLB(mu)
        return normHatE / np.sqrt(alp)



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
        return u_proj
        
    def RB_toNumpy(self):
        """Return the basis as a NumPy matrix
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
        print("[reducedBasis] Start greedy algorithm")
        self.N = 0
        S = []
        self.Z = []
        Delta = 1 + eps_tol

        self.DeltaMax = []

        self.computeOfflineErrorRhs()

        mu = mu_0

        # fig = go.Figure()
        
        while Delta > eps_tol and self.N < Nmax:

            S.append(mu)

            # self.computeOfflineReducedBasis(S)
            # self.computeOfflineError()

            if self.N == 0:
                self.computeOfflineReducedBasis(S)
                self.computeOfflineError()
            else:
                beta = self.model.computeBetaQm(mu)
                A = self.assembleA(beta[0][0])
                F = self.assembleF(beta[1][0][0])
                
                pc = self.ksp.getPC()
                pc.setType(self.PC_TYPE)

                self.reshist = {}

                self.ksp.setOperators(A)
                self.ksp.setConvergenceHistory()
                sol = F.duplicate()
                self.ksp.solve(F, sol)
                self.Z.append(sol)
                self.orthonormalizeZ()

                self.expandOffline()

                self.N += 1

                self.generateANq()
                self.generateLNp()
                self.generateFNp()

            mu_max = 0
            i_max = 0
            Delta_max = -float('inf')

            # l = []
            # m = []

            for i,mu_tmp in enumerate(tqdm(Dmu, desc=f"[reducedBasis] Greedy, step {self.N}", ascii=False, ncols=120)):
                beta = self.model.computeBetaQm(mu_tmp)
                ANmu = self.assembleAN(beta[0][0])
                uN,_ = self.getSolutions(mu_tmp)
                norm_uMu = np.sqrt( uN.T @ ANmu @ uN )

                precalc = {"betaA":beta[0][0], "betaF":beta[1][0][0], "uN":uN}

                D_en = self.computeEnergyBound(mu_tmp, precalc=precalc)
                Delta_tmp = D_en / norm_uMu

                # m.append(mu_tmp.parameterNamed('k_1'))
                # l.append(Delta_tmp)

                if Delta_tmp > Delta_max:
                    i_max = i
                    Delta_max = Delta_tmp
                    mu_max = mu_tmp

            Dmu.pop(i_max)
            Delta = Delta_max
            self.DeltaMax.append(Delta_max)
            mu = mu_max
            print("[reducedBasis] Greedy alg.", self.N, Delta, mu)

        
        if self.N == Nmax:
            print("[reducedBasis] Greedy alg : warning, max size reached")

    #     fig.update_layout(
    #         # title="Energy for ",
    #         xaxis_title=r"$\mu$",
    #         yaxis_title=r"$\Delta_N$",
    #         legend_title=r"$N$"
    #     )
    #     fig.write_image("energy.png", scale=2)
    #     fig.write_html("energy.html")
        print("[reducedBasis] End greedy algorithm")
        self.DeltaMax = np.array(self.DeltaMax)

        return S



    """
    Save and load results
    """
    def saveReducedBasis(self, path, force=False):
        """save the reduced basis in files

        Args:
            path (str): path of the directory whre data are to be saved
            force (bool, optional): Force saving, even if files are already present. Defaults to False.
        """
        if os.path.isdir(path) and not force:
            print(f"[reducedBasis] Directory {path} already exists. Rerun with force=True to force saving")
            return
        elif not os.path.isdir(path):
            os.mkdir(path)
        
        if not os.path.isdir(path+'/ANq'):
            os.mkdir(path+'/ANq')
        if not os.path.isdir(path+'/FNp'):
            os.mkdir(path+'/FNp')
        if not os.path.isdir(path+'/offline'):
            os.mkdir(path+'/offline')

        for q,Aq in enumerate(self.ANq):
            np.save(path+"/ANq/"+str(q), Aq)
        for p,Fp in enumerate(self.FNp):
            np.save(path+"/FNp/"+str(p), Fp)


        np.save(path+'/offline/SS', self.SS)
        np.save(path+'/offline/SL', self.SL)
        np.save(path+'/offline/LL', self.LL)

        if self.DeltaMax is not None:
            np.save(path+'/offline/DeltaMax', self.DeltaMax)

    def loadReducedBasis(path, model):
        if not os.path.isdir(path):
            print(f"[reducedBasis] Error : could not find {path}")
            return None
        if not os.path.isdir(path + '/ANq'):
            print(f"[reducedBasis] Error : could not find {path}/ANq")
            return None
        if not os.path.isdir(path + '/FNp'):
            print(f"[reducedBasis] Error : could not find {path}/FNp")
            return None
        if not os.path.isdir(path + '/offline'):
            print(f"[reducedBasis] Error : could not find {path}/offline")
            return None
        
        rbLoaded = reducedbasis(None, None, model, None, None)
        rbLoaded.Qa = len(os.listdir(path+'/ANq'))
        rbLoaded.Qf = len(os.listdir(path+'/FNp'))
        rbLoaded.ANq = []
        rbLoaded.FNp = []

        for q in range(rbLoaded.Qa):
            rbLoaded.ANq.append(np.load(path+'/ANq/'+str(q)+'.npy'))
        for p in range(rbLoaded.Qf):
            rbLoaded.FNp.append(np.load(path+'/FNp/'+str(p)+'.npy'))
        rbLoaded.SS = np.load(path+'/offline/SS.npy')
        rbLoaded.SL = np.load(path+'/offline/SL.npy')
        rbLoaded.LL = np.load(path+'/offline/LL.npy')

        rbLoaded.N = rbLoaded.ANq[0].shape[0]

        if os.path.isfile(path+'/offline/DeltaMax.npy'):
            rbLoaded.DeltaMax = np.load(path+'/offline/DeltaMax.npy')

        return rbLoaded
