from .reducedbasis import *
import scipy.linalg as sl
import warnings


def taille_dict(d):
    n = 0
    for k in d:
        n += len(d[k])
    return n

class reducedbasisTime(reducedbasis):

    def __init__(self, *, tf, K, **kwargs) -> None:
        """Initialize the object

        Args:
            tf (float) : final time
            K (int)    : number of time iterations
            **kwargs   : keyword arguments for the reducedbasis class
        """

        warnings.warn("This class is still in construction. Correction from the previous version are in progress")

        super().__init__(**kwargs)
        
        self.Qm = 0         # size of the decomposition of the mass matrix M

        self.tf = tf        # final time
        self.K  = K         # number of time iterations
        self.dt = tf / K    # time step

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

        # MN = np.einsum('r,rnm->nm', beta, self.MNr)
        return self.MNr[0]


    """Solve for a given parameter
    """

    def solveTime(self, mu, g, beta=None):
        """Solve the time-dependent equation for a given parameter and time-dependent function

        Args:
            mu (ParameterSpaceElement): parameter used
            g (function): right-hand side time-dependent function
            beta (list, optional) : coefficients of the decomposition, if they have already been computed


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
        u = np.zeros(self.N)    # initial solution TODO

        for k in range(1, self.K):
            sol = sl.lu_solve(matLu, g(k * self.dt) * self.dt * FNmu + MNmu @ u)
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
        """Compute the online error SQUARED, from a parameter mu, at time k

        Args:
            mu (ParameterSpaceElement): parameter used
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


    def computeOnlineError(self, mu, beta, g, computeEnergyNorm=False):
        """Compute online bound error, from a parameter mu

        Args:
            mu (ParameterSpaceElement): parameter used
            beta (list) : coefficients of the decomposition
            g (function): right-hand side time-dependent function
            computeEnergyNorm (bool): computes the energy normsuring the resolution (stroed in self.EnNorm).\
                Defaults to False

        Returns:
            float: Δ_N^k
        """
        betaA = beta[0][0]
        betaF = beta[1][0][0]
        betaM = beta[2][0]
        ANmu = self.assembleAN(betaA)
        FNmu = self.assembleFN(betaF)
        MNmu = self.assembleMN(betaM)
        precalc = {
            "betaA": betaA,
            "betaF": betaF,
            "betaM": betaM
        }

        mat = MNmu + self.dt * ANmu
        matLu = sl.lu_factor(mat)

        alpm1 = 1 / np.sqrt(self.alphaLB_(beta))

        aNorm = np.zeros(self.K)
        self.err[:] = 0
        self.DeltaN[:] = 0

        u = np.zeros(self.N)    # initial condition TODO
        if computeEnergyNorm:
            self.EnNorm[:] = 0
        
        for k in range(1, self.K+1):
            u_tmp = sl.lu_solve(matLu, g(k * self.dt) * self.dt * FNmu + MNmu @ u )
            precalc["uN"] = u_tmp
            precalc["uNm1"] = u
            self.err[k-1] = self.computeOnlineError_k(mu, u_tmp, u, k-1, precalc=precalc)
            u = np.array(u_tmp).flatten()

            self.DeltaN[k-1] = alpm1 * np.sqrt(self.dt * self.err[:k].sum())

            if computeEnergyNorm:
                aNorm[k-1] = u.T @ ANmu @ u
                self.EnNorm[k-1] = np.sqrt( u.T @ MNmu @ u + self.dt * aNorm.sum() )

        return self.DeltaN[-1]


    """
    Save and load functions
    """
    def saveReducedBasis(self, path, force=False, check=True):
        """Save the reduced basis in files

        Args:
            path (str): path of the directory whre data are to be saved
            force (bool, optional): Force saving, even if files are already present. Defaults to False.
            check (bool, optional): Check that the exported values are correct (only in sequential). Defaults to True
        """
        h5f, content = super().saveReducedBasis(path, force=force, notDoneYet=True)
        jsonPath = f"{os.getcwd()}/reducedbasis.json"

        if self.worldComm.isMasterRank():

            content['Qm'] = self.Qm
            content['tf'] = self.tf
            content['K'] = self.K

            f = open(jsonPath, 'w')
            json.dump(content, f, indent=4)
            f.close()

            h5f.create_dataset("MNr", data=self.MNr)

            h5f.close() 
            print("Done !")

            if check and feelpp.Environment.isSequential():
                print("[reducedbasis] Checking that the exported basis is correct...")
                self.checkSaved(jsonPath)
                print("[reducedbasis] Check is ok !")


    def loadReducedBasis(self, path, model):
        """Load reduced basis from json

        Args:
            path (str): path to the json description file
            model (toolboxmor): toolboxmor used to create the model
        """
        h5f, j = super().loadReducedBasis(path, model, notDoneYet=True)

        self.tf = j['tf']
        self.K = j['K']
        self.dt = self.tf / self.K

        self.MNr = h5f['MNr'][:]

        h5f.close()