"""reducedBasisPetsc :
    deals with reduced basis, using PETSc matrices
"""

import sys, slepc4py
import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc
slepc4py.init(sys.argv)
# import plotly.graph_objects as go
import sys, os
import h5py
import json
from tqdm import tqdm
from scipy.sparse.linalg import splu, spsolve

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



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

    def __init__(self, model) -> None:
        """Initialise the object

        Args:
            model (ToolboxMor_{2|3}D): model DEIM used for the decomposition
        """
        self.Qa = 0
        self.Qf = 0
        self.N = 0

        self.model = model  # TODO : load online model
        
        self.ANq    : np.ndarray   # tensor of shape (Qa, N, N)
        self.FNp    : np.ndarray   # tensor of shape (Qf, N)

        self.SS = np.zeros((self.Qf, self.Qf))            # SS[p,p_] = (Sp,sp_)X
        self.SL     : np.ndarray    # shape (Qf, Qa, N)     SL[p,n,q] = (Sp,Lnq)X
        self.LL     : np.ndarray    # shape (Qa, N, Qa, N)  LL[n,p,n_,p_] = (Lnq,Ln_p_)X


        ## For greedy memory
        self.DeltaMax = None
        if rank == 0:
            print("[reducedbasis] rb initialized")

        def alphaLB(mu):
            if rank == 0:
                print("[reducedbasis] WARNING : Lower bound not initialized yet")
            return 1

        def alphaLB_(beta):
            if rank == 0:
                print("[reducedbasis] WARNING : Lower bound not initialized yet")
            return 1

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_

        self.isInitilized = False


    def setMubar(self, mubar):
        """Set mubar form a given parameter

        Args:
            mubar (ParameterSpaceElement) : parameter bar{µ}
        """
        self.mubar = mubar
        beta = self.model.computeBetaQm(self.mubar)
        self.betaA_bar = beta[0]
        self.betaF_bar = beta[1]


    """
    Online computation
    """
    def assembleAN(self, beta, size=None):
        """Assemble the reduced matrix from a given vector of parameters
           AN = sum_q beta[q]*ANq[q]

        Args:
            beta (list): list of parameters betaA
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            np.ndarray: assembled matrix
        """
        assert( len(beta) == self.Qa ), f"Number of param ({len(beta)}) should be {self.Qa}"

        AN = np.einsum('q,qnm->nm', beta, self.ANq, optimize=True)
        # AN = self.ANq[0].copy()
        # AN.fill(0)
        # for q in range(0, self.Qa):
        #     AN += self.ANq[q] * beta[q]

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
            np.ndarray: assembled vector
        """
        assert( len(beta) == self.Qf ), f"Number of param ({len(beta)}) should be {self.Qf}"

        FN = np.einsum('p,pn->n', beta, self.FNp, optimize=True)
        # FN = self.FNp[0].copy()
        # FN.fill(0)
        # for q in range(0, self.Qf):
        #     FN += self.FNp[q] * beta[q]

        if size is None:
            return FN
        else:
            return FN[:size]


    def getSolutions(self, mu, beta = None, size=None):
        """Return solution uN and output sN

        Args:
            mu (ParameterSpaceElement) : parameter
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            tuple (np.ndarray,float) : (uN, sN)
        """
        if beta is None:
            beta = self.model.computeBetaQm(mu)
        A_mu = self.assembleAN(beta[0][0], size=size)
        F_mu = self.assembleFN(beta[1][0][0], size=size)

        sol = np.linalg.solve(A_mu, F_mu)

        return sol, sol @ F_mu



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
        s2 = np.einsum('q,p,n,qpn', betaA, betaF, uN, self.SL, optimize=True)
        s3 = np.einsum('q,r,n,m,qnrm', betaA, betaA, uN, uN, self.LL, optimize=True)

        # s1 = 0
        # for p in range(self.Qf):
        #     for p_ in range(self.Qf):
        #         s1 += betaA[p] * betaA[p_] * self.SS[p,p_]
        # s2 = 0
        # for p in range(self.Qf):
        #     for q in range(self.Qa):
        #         for n in range(self.N):
        #             s2 += betaF[p] * betaA[q] * uN[n] * self.SL[q,p,n]
        #
        # s3 = 0
        # for q in range(self.Qa):
        #     for n in range(self.N):
        #         for q_ in range(self.Qa):
        #             for n_ in range(self.N):
        #                 s3 += betaA[q] * betaA[q_] * uN[n] * uN[n_] * self.LL[q,n,q_,n_]

        return s1 + 2*s2 + s3



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
        if precalc is None:
            alp = self.alphaLB(mu)
        else:
            alp = self.alphaLB_(precalc["betaA"])
        return normHatE / np.sqrt(alp)




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

        print(f"[reducedbasis] saving reduced basis to {os.getcwd()}/reducedbasis.json...", end=" ")

        content = {'Qa': self.Qa, 'Qf': self.Qf, 'N': self.N, "path": path+"/reducedbasis.h5"}
        os.system('pwd')
        h5f = h5py.File(path+"/reducedbasis.h5", "w")
        f = open('reducedbasis.json', 'w')

        dict_mubar = {}
        for n in self.mubar.parameterNames(): dict_mubar[n] = self.mubar.parameterNamed(n)
        content["mubar"] = dict_mubar

        h5f.create_dataset("ANq", data=self.ANq)
        h5f.create_dataset("FNp", data=self.FNp)
        # for q, Aq in enumerate(self.ANq):
        #     h5f.create_dataset(f"AN{q}", data=Aq)
        # for p, Fp in enumerate(self.FNp):
        #     h5f.create_dataset(f"FN{p}", data=Fp)
        h5f.create_dataset(f"SS", data=self.SS)
        h5f.create_dataset(f"SL", data=self.SL)
        h5f.create_dataset(f"LL", data=self.LL)

        if self.DeltaMax is not None:
            h5f.create_dataset("DeltaMax", data=self.DeltaMax)
        else:
            h5f.create_dataset("DeltaMax", data=np.array([]))

        h5f.create_dataset("alphaMubar", data=np.array([self.alphaMubar]))


        json.dump(content, f, indent = 4)
        h5f.close()
        f.close()
        print("Done !")
        return os.getcwd() + "/reducedbasis.json"


    def loadReducedBasis(self, path, model):
        """Load reduced basis from json

        Args:
            path (str): path to the json description file
            model (toolboxmor): toolboxmor used to create the model
        """
        f = open(path, "r")
        j = json.load(f)
        f.close()

        self.model = model

        self.N = j['N']
        self.Qa = j['Qa']
        self.Qf = j['Qf']
        mubar = model.parameterSpace().element()
        mubar.setParameters(j["mubar"])
        self.setMubar(mubar)

        h5f = h5py.File(j['path'], "r")

        self.ANq = h5f["ANq"][:]
        self.FNp = h5f["FNp"][:]

        self.SS = h5f["SS"][:]
        self.SL = h5f["SL"][:]
        self.LL = h5f["LL"][:]

        try:
            assert self.ANq.shape == (self.Qa, self.N, self.N), f"Wrong shape for ANq (excepted {(self.Qa, self.N, self.N)}, got {self.ANq.shape}"
            assert self.FNp.shape == (self.Qf, self.N), f"Wrong shape for ANq (excepted {(self.Qf, self.N)}, got {self.FNp.shape}"
            assert self.SS.shape == (self.Qf, self.Qf), f"Wrong shape for SS (excepted {(self.Qf, self.Qf)}, got {self.SS.shape}"
            assert self.SL.shape == (self.Qa, self.Qf, self.N), f"Wrong shape for SL (excepted {(self.Qa, self.Qf, self.N)}, got {self.SL.shape}"
            assert self.LL.shape == (self.Qa, self.N, self.Qa, self.N), f"Wrong shape for LL (excepted {(self.Qa, self.N, self.Qa, self.N)}, got {self.LL.shape}"
        except AssertionError as e:
            print(f"[reduced basis] Something went wrong when loading {path} : {e}")


        tmpDeltaMax = h5f["DeltaMax"][:]
        self.DeltaMax = None if tmpDeltaMax.shape == (0,) else tmpDeltaMax

        self.alphaMubar = h5f["alphaMubar"][0]
        betaA_bar_np = np.array(self.betaA_bar)
        def alphaLB(mu):
            # From a parameter
            betaMu = self.model.computeBetaQm(self.mu)[0][0]
            return self.alphaMubar * np.min( betaMu / betaA_bar_np )

        def alphaLB_(betaA):
            # From a decomposition
            return self.alphaMubar * np.min( betaA / betaA_bar_np )

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_

        print(f"[reduced basis] Basis loaded from {path}")
        self.setInitialized = True



class reducedbasisOffline(reducedbasis):
    """ Reduced basis for stationnary problem
    """

    
    def __init__(self, Aq, Fq, model, mubar):

        super().__init__(model)

        self.setMubar(mubar)

        self.Aq = Aq
        self.Fq = Fq

        self.Qa = len(Aq)
        self.Qf = len(Fq)


        self.SS = np.zeros((self.Qf, self.Qf))

        self.NN = Aq[0].size[0]

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
            if rank == 0:
                print("[petsc4py] Iteration {} Residual norm: {}".format(its,rnorm))
        
        def monitor_nverbose(ksp, its, rnorm):
            self.reshist[its] = rnorm
        
        self.monitorFunc = [monitor_nverbose, monitor]

        self.ksp.setMonitor(monitor_nverbose)
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
            # From a parameter
            betaMu = self.model.computeBetaQm(self.mu)[0][0]
            return alphaMubar * np.min( betaMu / betaA_bar_np )

        def alphaLB_(betaA):
            # From a decomposition
            return alphaMubar * np.min( betaA / betaA_bar_np )

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_



    
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
        elif rank == 0:
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

        for i,mu in enumerate(tqdm(mus, desc=f"[reducedBasis] Offline generation of the basis", ascii=False, ncols=120)):
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
        # self.ANq = []
        # for q in range(self.Qa):
        #     self.ANq.append(np.zeros((self.N,self.N)))
        self.ANq = np.zeros((self.Qa, self.N, self.N))
        for i,u in enumerate(self.Z):
            for j,v in enumerate(self.Z):
                for q in range(self.Qa):
                    self.ANq[q,i,j] = v.dot( self.Aq[q] * u )
                    # print(i, j, self.ANq[q][i,j]) # uT.A.u
    
    def expandANq(self):
        """Expand the reduced matrices ANq
        """
        self.ANq = np.concatenate( (self.ANq, np.zeros((self.Qa, self.N-1, 1))), axis=2)    # N has already been increased
        self.ANq = np.concatenate( (self.ANq, np.zeros((self.Qa, 1, self.N))), axis=1)
        for q in range(self.Qa):
            for i,u in enumerate(self.Z):
                self.ANq[q, -1, i] = u.dot( self.Aq[q] * self.Z[-1] )
                self.ANq[q, i, -1] = self.Z[-1].dot( self.Aq[q] * u )
    
    def generateLNp(self):
        """Generate the reduced vectors LNq
        """
        self.LNp = np.zeros((self.Qf, self.N))
        for i,u in enumerate(self.Z):
            for q in range(self.Qf):
                self.LNp[q,i] = self.Fq[q].dot(u)
    
    def generateFNp(self):
        """Generate the reduced vectors FNp
        """
        self.FNp = np.zeros((self.Qf, self.N))
        for i,u in enumerate(self.Z):
            for q in range(self.Qf):
                self.FNp[q,i] = self.Fq[q].dot(u)
    
    def expandFNp(self):
        self.FNp = np.concatenate( (self.FNp, np.zeros((self.Qf,1))), axis=1)
        for q in range(self.Qf):
            self.FNp[q, -1] = self.Fq[q].dot(self.Z[-1])

    def expandLNp(self):
        self.LNp = np.concatenate( (self.LNp, np.zeros((self.Qf,1))), axis=1)
        for q in range(self.Qf):
            self.LNp[q, -1] = self.Fq[q].dot(self.Z[-1])



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
        if rank == 0:
            print("[reducedBasis] Start greedy algorithm")
        self.N = 0
        S = []
        self.Z = []
        Delta = 1 + eps_tol

        self.DeltaMax = []

        self.computeOfflineErrorRhs()

        mu = mu_0

        betas = {}
        for mu in Dmu:
            betas[mu] = self.model.computeBetaQm(mu)

        # fig = go.Figure()
        
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
                beta = betas[mu_tmp]
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

            if rank == 0:
                print(f"[reducedBasis] Greedy algo, N={self.N}, Δ={Delta:e} (tol={eps_tol:e})".ljust(64), f"µ={mu}")
        if rank == 0 and self.N == Nmax:
            print("[reducedBasis] Greedy algo : warning, max size reached")

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

