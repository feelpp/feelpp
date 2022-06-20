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
# import time

import feelpp




def convertToPetscMat(Aq):
    """convert a list of feelpp._alg.MatrixPetscDouble to a list of petsc4py.PETSc.Mat
    """
    AqP = []
    for A in Aq:
        AqP.append(A.to_petsc().mat())
        if not AqP[-1].assembled:
            AqP[-1].assemble()
    return AqP

def convertToPetscVec(Fq):
    """convert a list of feelpp._alg.VectorPetscDouble to a list of petsc4py.PETSc.Vec
    """
    FqP = []
    for F in Fq:
        FqP.append(F.to_petsc().vec())
    return FqP


class reducedbasis():

    def __init__(self, model, worldComm=None) -> None:
        """Initialise the object

        Args:
            model (ToolboxMor_{2|3}D): model DEIM used for the decomposition
        """
        self.Qa = 0                 # size of the decomposition of A
        self.Qf = 0                 # size of the decomposition of F
        self.QLk = []               # sizes of the decompositions of Lk
        self.N_output = 0           # number of outputs 
        self.output_names : list    # names of the outputs
        
        self.N = 0      # size of the rediced basis

        if worldComm == None:
            self.worldComm = feelpp.Environment.worldCommPtr()
        else:
            self.worldComm = worldComm

        self.model = model  # TODO : load online model

        self.ANq    : np.ndarray   # tensor of shape (Qa, N, N)
        self.FNp    : np.ndarray   # tensor of shape (Qf, N)
        self.LkNp   : list         # list of len n_output of tensors of shape (QLk, N)


        self.SS = np.zeros((self.Qf, self.Qf))            # SS[p,p_] = (Sp,sp_)X
        self.SL     : np.ndarray    # shape (Qf, Qa, N)     SL[p,n,q] = (Sp,Lnq)X
        self.LL     : np.ndarray    # shape (Qa, N, Qa, N)  LL[n,p,n_,p_] = (Lnq,Ln_p_)X


        ## For greedy memory
        self.DeltaMax = None

        def alphaLB(mu):
            if self.worldComm.isMasterRank():
                print("[reducedbasis] WARNING : Lower bound not initialized yet")
            return 1

        def alphaLB_(beta):
            if self.worldComm.isMasterRank():
                print("[reducedbasis] WARNING : Lower bound not initialized yet")
            return 1

        def gammaUB(mu):
            if self.worldComm.isMasterRank():
                print("[reducedbasis] WARNING : Upper bound not initialized yet")
            return 1

        def gammaUB_(beta):
            if self.worldComm.isMasterRank():
                print("[reducedbasis] WARNING : Upper bound not initialized yet")
            return 1

        self.alphaMubar = 1
        self.gammaMubar = 1

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_
        self.gammaUB = gammaUB
        self.gammaUB_ = gammaUB_

        self.isInitilized = False
        self.mubar = None

        if self.worldComm.isMasterRank():
            print("[reducedbasis] Online rb initialized")


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

    def assembleLkN(self, k, beta, size=None):
        """Assemble the rhs from a given vector of parameters
           Lk = sum_q beta[q]*Lkq[q]

        Args:
            k (int): number of the i-th output
            beta (list): list of parameters betaF

        Returns:
            PETSc.Vec: assemble vector
        """

        assert( len(beta) == self.QLk[k]), f"Number of param ({len(beta)}) should be {self.QLk[k]}"

        LkN = np.einsum('p,pn->n', beta, self.LkNp[k])

        if size is None:
            return LkN
        else:
            return LkN[:size]




    def getSolutions(self, mu, beta=None, k=-1, size=None):
        """Return solution uN and output sN

        Args:
            mu (ParameterSpaceElement): parameter
            beta (list, optional) : coefficients of the decomposition, if they have already been computed
            k (int, optional) : index of the output to be computed, if -1 the compliant output is computed
            size (int, optional): size of the sub-basis wanted. Defaults to None, meaning the whole basis is used.

        Returns:
            tuple (np.ndarray,float) : (uN, sN)
        """
        if beta is None:
            beta = self.model.computeBetaQm(mu)
        A_mu = self.assembleAN(beta[0][0], size=size)
        F_mu = self.assembleFN(beta[1][0][0], size=size)

        sol = np.linalg.solve(A_mu, F_mu)

        if k == -1:
            l = F_mu
        else:
            if 0 <= k and k < self.N_output:
                l = self.assembleLkN(k, beta[1][k+1][0])
            else:
                raise ValueError(f"Output {k} not valid")

        return sol, sol @ l


    def getOutputName(self, k):
        """Get the name of the k-th output

        Args:
            k (int): index of the output

        Returns:
            str: name of the output
        """
        if k == -1:
            return "Compliant"
        else:
            return self.output_names[k]




    def computeOnlineError(self, mu, precalc=None, size=None):
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
            A_mu = self.assembleAN(betaA, size=size)
            F_mu = self.assembleFN(betaF, size=size)

            uN = np.linalg.solve(A_mu, F_mu)
        else:
            betaA = precalc["betaA"]
            betaF = precalc["betaF"]
            uN = precalc["uN"]
        s = [size, self.N][size is None]

        s1 = betaF @ self.SS @ betaF    # beta_p*beta_p'*(Sp,Sp')
        s2 = np.einsum('q,p,n,qpn', betaA, betaF, uN, self.SL[:,:,:s], optimize=True)
        s3 = np.einsum('q,r,n,m,qnrm', betaA, betaA, uN, uN, self.LL[:,:s,:,:s], optimize=True)

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



    def computeEnergyBound(self, mu, precalc=None, size=None):
        """computes ernergy bound

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function
            size (int, optional): size of the sub-basis considered

        Returns:
            float: energy error bound
        """
        normHatE = self.computeOnlineError(mu, precalc=precalc, size=size)
        if precalc is None:
            alp = self.alphaLB(mu)
        else:
            alp = self.alphaLB_(precalc["betaA"])
        return normHatE / np.sqrt(alp)

    def computeEffectivity(self, mu, precalc=None, size=None):
        """computes effectivity = Delta_N / ||e(mu)|_X

        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function
            size (int, optional): size of the sub-basis considered

        Returns:
            float: effectivity
        """
        normHatE = self.computeOnlineError(mu, precalc=precalc, size=size)
        if precalc is None:
            alp = self.alphaLB(mu)
        else:
            alp = self.alphaLB_(precalc["betaA"])
        return normHatE / alp


    def UB_LB(self, mu, precalc=None):
        """computed the theorical bound for effectivity (1 ⩽ η(µ) ⩽ γ(µ)/α(µ) )
                                                                    ^^^^^^^^^^
        Args:
            mu (ParameterSpaceElement): parameter
            precalc (dict, optional): Dict containing the values of betaA, betaF and uN if these values have already been calculated. Defaults to None.\
                If None is given, the quantities are calculated in the function

        Return:
            float: gamma(mu)/alpha(mu)
        """
        if precalc is None:
            betaMu = self.model.computeBetaQm(mu)[0][0]
        else:
            betaMu = precalc["betaA"]
        gamma = self.gammaUB_(betaMu)
        alpha = self.alphaLB_(betaMu)
        return gamma/alpha


    """
    Save and load results
    """
    def saveReducedBasis(self, path, force=False, check=True):
        """save the reduced basis in files

        Args:
            path (str): path of the directory whre data are to be saved
            force (bool, optional): Force saving, even if files are already present. Defaults to False.
            check (bool, optional): Check that the exported values are correct (only in sequential). Defaults to True
        """
        if os.path.isdir(path) and not force:
            if self.worldComm.isMasterRank():
                print(f"[reducedBasis] Directory {path} already exists. Rerun with force=True to force saving")
            return
        elif not os.path.isdir(path):
            if self.worldComm.isMasterRank():
                os.makedirs(path, exist_ok=True)

        if self.worldComm.isMasterRank():
            print(f"[reducedbasis] saving reduced basis to {os.getcwd()}/reducedbasis.json ...", end=" ")

        jsonPath = f"{os.getcwd()}/reducedbasis.json"

        if self.worldComm.isMasterRank():
            h5f = h5py.File(path+"/reducedbasis.h5", "w")
            content = {'Qa': self.Qa, 'Qf': self.Qf, 'N_output': self.N_output, "QLk": self.QLk, "output_names": self.output_names, 'N': self.N, "path": path+"/reducedbasis.h5"}
            dict_mubar = {}
            for n in self.mubar.parameterNames(): dict_mubar[n] = self.mubar.parameterNamed(n)
            content["mubar"] = dict_mubar

            f = open('reducedbasis.json', 'w')
            json.dump(content, f, indent = 4)
            f.close()


            h5f.create_dataset("ANq", data=self.ANq)
            h5f.create_dataset("FNp", data=self.FNp)
            for k in range(self.N_output):
                h5f.create_dataset(f"L{k}Np", data=self.LkNp[k])

            h5f.create_dataset("SS", data=self.SS)
            h5f.create_dataset("SL", data=self.SL)
            h5f.create_dataset("LL", data=self.LL)

            if self.DeltaMax is not None:
                h5f.create_dataset("DeltaMax", data=self.DeltaMax)
            else:
                h5f.create_dataset("DeltaMax", data=np.array([]))

            h5f.create_dataset("alphaMubar", data=np.array([self.alphaMubar]))
            h5f.create_dataset("gammaMubar", data=np.array([self.gammaMubar]))

            h5f.close()

            print("Done !")

            if check and feelpp.Environment.isSequential():
                print("[reducedbasis] Checking that the exported basis is correct...")
                self.checkSaved(jsonPath)
                print("[reducedbasis] Check is ok !")

        return f"{os.getcwd()}/reducedbasis.json"


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
        self.N_output = j['N_output']
        self.QLk = j['QLk']
        self.output_names = j['output_names']
        mubar = model.parameterSpace().element()
        mubar.setParameters(j["mubar"])
        self.setMubar(mubar)

        h5f = h5py.File(j['path'], "r")

        self.ANq = h5f["ANq"][:]
        self.FNp = h5f["FNp"][:]
        self.LkNp = []
        for k in range(self.N_output):
            Lkp = h5f[f"L{k}Np"][:]
            self.LkNp.append(Lkp)

        self.SS = h5f["SS"][:]
        self.SL = h5f["SL"][:]
        self.LL = h5f["LL"][:]

        try:
            assert self.ANq.shape == (self.Qa, self.N, self.N), f"Wrong shape for ANq (excepted {(self.Qa, self.N, self.N)}, got {self.ANq.shape}"
            assert self.FNp.shape == (self.Qf, self.N), f"Wrong shape for FNp (excepted {(self.Qf, self.N)}, got {self.FNp.shape}"
            for k in range(self.N_output):
                assert self.LkNp[k].shape == (self.QLk[k], self.N), f"Wrong shape for L{k}Np (excepted {(self.QLk[k], self.N)}, got {self.LkNp[k].shape}"

            assert self.SS.shape == (self.Qf, self.Qf), f"Wrong shape for SS (excepted {(self.Qf, self.Qf)}, got {self.SS.shape}"
            assert self.SL.shape == (self.Qa, self.Qf, self.N), f"Wrong shape for SL (excepted {(self.Qa, self.Qf, self.N)}, got {self.SL.shape}"
            assert self.LL.shape == (self.Qa, self.N, self.Qa, self.N), f"Wrong shape for LL (excepted {(self.Qa, self.N, self.Qa, self.N)}, got {self.LL.shape}"
        except AssertionError as e:
            print(f"[reduced basis] [{self.worldComm.localRank()}] Something went wrong when loading {path} : {e}")


        tmpDeltaMax = h5f["DeltaMax"][:]
        self.DeltaMax = None if tmpDeltaMax.shape == (0,) else tmpDeltaMax

        self.alphaMubar = h5f["alphaMubar"][0]
        self.gammaMubar = h5f["gammaMubar"][0]
        betaA_bar_np = np.array(self.betaA_bar)

        def alphaLB(mu):
            # From a parameter
            betaMu = self.model.computeBetaQm(mu)[0][0]
            return self.alphaMubar * np.min( betaMu / betaA_bar_np )

        def alphaLB_(betaA):
            # From a decomposition
            return self.alphaMubar * np.min( betaA / betaA_bar_np )

        def gammaUB(mu):
            # From a parameter
            betaMu = self.model.computeBetaQm(mu)[0][0]
            return self.gammaMubar * np.max( betaMu / betaA_bar_np )

        def gammaUB_(betaA):
            # From a decomposition
            return self.gammaMubar * np.max( betaA / betaA_bar_np )

        self.alphaLB = alphaLB
        self.alphaLB_ = alphaLB_
        self.gammaUB = gammaUB
        self.gammaUB_ = gammaUB_

        if self.worldComm.isMasterRank():
            print(f"[reduced basis] Basis loaded from {path}")
        self.setInitialized = True


    def checkSaved(self, path):
        """Checks that the basis has been exported correctly

        Args:
            path (string): path to the reducedbasis.json file
        """

        rbLoaded = reducedbasis(None)
        rbLoaded.loadReducedBasis(path, self.model)

        for n in self.mubar.parameterNames():
            assert self.mubar.parameterNamed(n) == rbLoaded.mubar.parameterNamed(n)

        assert( self.Qa == rbLoaded.Qa )
        assert( self.Qf == rbLoaded.Qf )
        assert( self.N == rbLoaded.N )

        for q in range(rbLoaded.Qa):
            assert( (self.ANq[q] == rbLoaded.ANq[q]).all() ), f"q = {q}"
        for p in range(rbLoaded.Qf):
            assert( (self.FNp[p] == rbLoaded.FNp[p]).all() ), f"p = {p}"

        assert( (self.SS == rbLoaded.SS).all() ), "SS"
        assert( (self.SL == rbLoaded.SL).all() ), "SL"
        assert( (self.LL == rbLoaded.LL).all() ), "LL"

        if self.DeltaMax is None:
            assert rbLoaded.DeltaMax is None, "DeltaMax"
        else:
            assert( (self.DeltaMax == rbLoaded.DeltaMax).all() ), "DeltaMax"

        assert self.alphaMubar == rbLoaded.alphaMubar, "alphaMubar"
        assert self.gammaMubar == rbLoaded.gammaMubar, "gammaMubar"


    def getSize(self, threshold):
        """Get size of the basis for a desired threshold error

        Args:
            threshold (float): thresohld

        Returns:
            int: size of the basis
        """
        if self.DeltaMax is None:
            if self.worldComm.isMasterRank():
                print("[reducedbasis] Maximal error has not been calculated yet")
            return self.N
        ind = 0
        while self.DeltaMax[ind] > threshold:
            ind += 1
        return ind

