from .reducedbasis import *
import scipy.linalg as sl
import warnings


def taille_dict(d):
    n = 0
    for k in d:
        n += len(d[k])
    return n

class reducedbasisTime(reducedbasis):

    def __init__(self, Aq, Fp, model, mubar, alphaLB, Mr, tf, K) -> None:
        """Initialise the object

        Args:
            `Aq` (list of PETSc.Mat)       : matrices Aq given by the affine decomposition
            `Fp` (list of PETSc.Vec)       : vectors Fq given by the decomposition of right-hand side
            `model` (ToolboxMor_{2|3}D)    : model DEIM used for the decomposition
            `mubar` (ParameterSpaceElement): parameter mu_bar for the enrgy norm
            `alphaLB` (func)               : function mu ↦ alphaLB(mu)
            `Mr` (list of PETSc.Mat)       : matrices Mq of mass for the scalar product
            `tf` (float)                   : final time
            `K` (int)                      : number of time iterations
        """

        warnings("This class is still in construction. Correction from the previous version are in progress")

        super().__init__(model)
        self.Qm = len(Mr)
        # self.MN : spsp.csc_matrix

        self.tf = tf
        self.K  = K
        self.dt = tf / K

        self.MNr : np.ndarray # tensor of chape (Qm, N, N)


        self.FF = np.zeros((K, self.Qf, self.Qf)) #       FF[k,p,p_]    = (Fkp,Fkp_)X
        self.FM : np.ndarray    # shape (K, Qf, Qm, N)    FM[k,p,r,n]   = (Fkp, Mrn)X
        self.FL : np.ndarray    # shape (K, Qf, Qa, N)    FL[k,p,q,n]   = (Lnq, Fkp)X
        self.ML : np.ndarray    # shape (Qm, N, Qa, N)    ML[r,n,q,n_]  = (Lnq, Mrn_)X
        self.MM : np.ndarray    # shape (Qm, N, Qm, N)    MM[q,n,q_,n_] = (Mrn, Mr_n_)X

        # quantities depending on mu, but stored to avoid re-computation
        self.err = np.zeros(K)          # online error
        self.DeltaN = np.zeros(K)       # error bound using offline / online computation
        self.EnNorm = np.zeros(K)       # energy norm, filled as the resolution goes



    def assembleMN(self, beta):
        """Assemble the matrix from a given vector of parameters
           M = sum_r beta[r]*Mr[r]

        Args:
            beta (list): list of parameters betaM

        Returns:
            PETSc.Mat: assembled matrix
        """
        assert( len(beta) == self.Qm ), f"Number of param ({len(beta)}) should be {self.Qm}"

        MN = np.einsum('r,rnm->nm', beta, self.MNr)
        return MN


    """Solve for a given parameter
    """

    def solveTime(self, mu, g, beta=None):
        """Solves the time-dependant equation for a given parameter and time-dependant function

        Args:
            mu (ParameterSpaceElement): parameter
            g (np.ndarray): g function of the right-hand side (g[k]=g(k*dt))

        Returns:
            np.ndarray: solution uN of the equation at time t = K*dt
        """
        if beta is None:
            beta = self.model.computeBetaQm(mu)
        ANmu = self.assembleAN(beta[0][0])
        FNmu = self.assembleFN(beta[1][0][0])
        MNmu = np.array(self.assembleMN(beta[2][0]))

        mat = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(mat)
        u = np.zeros(self.N)    # initial solution

        for k in range(1, self.K):
            sol = sl.lu_solve(matLu, g[k] * self.dt * FNmu + MNmu @ u)
            u = sol.copy()
        return u



    """
    Error handling
    """
    def computeDirectError(self, mu, precalc=None):
        # TODO
        warnings.warn("The function computeDirectError has not yet been implemented")
        return super().computeDirectError(mu, precalc=precalc)



    def computeOnlineError_k(self, mu, uN, uNm1, k, precalc=None):
        """Computes the online error SQUARED, for a parameter mu

        Args:
            mu (ParameterSpaceElement): parameter
            uk (np.ndarray): solution at time k (of size N)
            ukm1 (np.ndarray): solution at time k-1
            k (int): index of current temps

        Returns:
            float: ||hat{e}^k||^2
        """
        if precalc is None:
            beta_ = self.model.computeBetaQm(mu)
            betaA = beta_[0][0]
            betaF = beta_[1][0][0]
            betaM = beta_[2]
            A_mu = self.assembleAN(betaA)
            F_mu = self.assembleFN(betaF)
            M_mu = self.assembleMN(betaM)

            uN = np.linalg.solve(A_mu, F_mu) # ?
        else:
            betaA = precalc["betaA"]
            betaF = precalc["betaF"]
            betaM = precalc["betaM"]
            uN = precalc["uN"]
            uNm1 = precalc["uNm1"]

        diff = ((uN.flatten() - uNm1.flatten()) / self.dt).T

        s1 = betaF @ self.FF[k] @ betaF
        s2 = np.einsum('p,r,n,prn', betaF, betaM, diff, self.FM[k])
        s3 = np.einsum('p,q,n,pqn', betaF, betaA, uN, self.FL[k])
        s4 = np.einsum('r,p,n,m,rnpm', betaM, betaA, diff, uN, self.ML)
        s5 = np.einsum('r,s,n,m,rnsm', betaM, betaM, diff, diff, self.MM)
        s6 = np.einsum('q,r,n,m,qnrm', betaA, betaA, uN, uN, self.LL)


        return s1 + 2 * (s2 + s3 + s4) + s5 + s6


    def computeOnlineError(self, mu, g, computeEnergyNorm=False):
        """Computes online bound error

        Args:
            mu (ParameterSpaceElement): parameter
            g (np.ndarray): right-hand side time-dependant function
            computeEnergyNorm (bool): computes the energy normsuring the resolution (stroed in self.EnNorm). Dafault to False

        Returns:
            float: Δ_N^k
        """
        beta_ = self.model.computeBetaQm(mu)
        betaA = beta_[0][0]
        betaF = beta_[1][0][0]
        betaM = beta_[2][0]
        ANmu = self.assembleAN(betaA)
        FNmu = self.assembleFN(betaF)
        MNmu = self.assembleMN(betaM)
        precalc = {
            "betaA": betaA,
            "betaF": betaF,
            "betaM": betaM
        }

        matN = splu(MNmu + self.dt * spsp.csc_matrix(ANmu))
        alpm1 = 1 / np.sqrt(self.alphaLB(mu))

        aNorm = np.zeros(self.K)
        self.err[:] = 0
        self.DeltaN[:] = 0

        u = np.zeros(self.N)    # initial condition
        if computeEnergyNorm:
            self.EnNorm[:] = 0
        
        for k in tqdm(range(1, self.K+1)):
            u_tmp = matN.solve( (g[k] * self.dt * FNmu + MNmu @ u).T )#.reshape(self.N)
            precalc["uN"] = u_tmp
            precalc["uNm1"] = u
            self.err[k-1] = self.computeOnlineError_k(mu, u_tmp, u, k-1, precalc=precalc)
            u = np.array(u_tmp).flatten()

            self.DeltaN[k-1] = alpm1 * np.sqrt(self.dt * self.err[:k].sum())

            if computeEnergyNorm:
                aNorm[k-1] = u.T @ ANmu @ u
                self.EnNorm[k-1] = np.sqrt( u.T @ MNmu @ u + self.dt * aNorm.sum() )

        return self.DeltaN[-1]

