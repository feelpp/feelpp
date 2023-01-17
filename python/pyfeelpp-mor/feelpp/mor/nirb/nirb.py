# -*- coding: utf-8 -*-
## Non intrusive reduced basis method with Feel++
## Thomas Saigre, Ali Elarif
## 09/2022


import feelpp
from .utils import *
import feelpp.operators as FppOp
# from NIRBinitCase import *
import numpy as np
import feelpp.toolboxes.heat as heat
import feelpp.toolboxes.fluid as fluid
import feelpp.interpolation as fi
from tqdm import tqdm
from random import choices
import math 
import pathlib


import os

from mpi4py import MPI

class ToolboxModel():
    """
        class containing the common tools associated with the toolbox used 
                    in the online and offline parts
    """
    def __init__(self, dim, H, h, toolboxType, model_path, finemesh_path, coarsemesh_path=None, order=1, **kwargs) -> None:
        """Initialize the toolbox model class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            finemesh_path (str): path to fine mesh file (if geo file, this will be the same as coarse mesh file)
            coarsemesh_path(str): path to coarse mesh file. Defaults to None. 
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """
        self.dimension = dim
        assert self.dimension in [2,3]
        self.H = H
        self.h = h if h != "H**2:H" else self.H**2
        self.order = order

        self.toolboxType = toolboxType
        assert self.toolboxType in ["heat", "fluid"], "toolboxType must be 'heat' or 'fluid'"
        self.model_path = feelpp.Environment.expand(model_path)
        self.finemesh_path = feelpp.Environment.expand(finemesh_path)

        self.worldcomm = feelpp.Environment.worldComm()

        if pathlib.Path(self.finemesh_path).suffix==".geo" :
            self.coarsemesh_path = self.finemesh_path
        else :
            self.coarsemesh_path = coarsemesh_path 

        assert self.coarsemesh_path != None, f"Set coarse mesh path"

        self.tbCoarse  = None
        self.tbFine    = None

        self.Xh = None          # fine function space

        self.outdir = None      # output directory

        self.initModel()
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Initialization done")


    def initModel(self):
        """Initialize the model
        """
        self.model = feelpp.readJson(self.model_path)
        self.tbFine = self.setToolbox(self.h)
        self.Xh = feelpp.functionSpace(mesh=self.tbFine.mesh(), order=self.order)
        self.Dmu = loadParameterSpace(self.model_path)
        self.Ndofs = self.getFieldSpace().nDof()

        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Number of nodes on the fine mesh : {self.tbFine.mesh().numGlobalPoints()}")

    def getFieldSpace(self, coarse=False):
        """Get the field space

        Args:
            coarse (bool, optional): get the coarse space. Defaults to False.

        Returns:
            feelpp._discr.*: field space
        """
        if not coarse:
            if self.toolboxType == "heat":
                return self.tbFine.spaceTemperature()
            elif self.toolboxType == "fluid":
                return self.tbFine.spaceVelocity()
            else:
                raise ValueError("Unknown toolbox")
        else:
            if self.toolboxType == "heat":
                return self.tbCoarse.spaceTemperature()
            elif self.toolboxType == "fluid":
                return self.tbCoarse.spaceVelocity()
            else:
                raise ValueError("Unknown toolbox")


    def initCoarseToolbox(self):
        """Initialize the rectification problem
        """
        self.tbCoarse = self.setToolbox(self.H,mesh_path=self.coarsemesh_path)
        self.XH = feelpp.functionSpace(mesh=self.tbCoarse.mesh(), order=self.order)
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Number of nodes on the coarse mesh : {self.tbCoarse.mesh().numGlobalPoints()}")


    def setToolbox(self, hsize, mesh_path=None):
        """Set up the toolbox object for the given model and mesh

        Args:
            hsize (float): size of the mesh

        Returns:
            Toolbox: toolbox object
        """
        if mesh_path==None: mesh_path=self.finemesh_path 
        # load meshes
        mesh_ = feelpp.mesh(dim=self.dimension, realdim=self.dimension)
        mesh = feelpp.load(mesh_, mesh_path, hsize)

        # set mesh and model properties
        if self.toolboxType == "heat":
            tb = heat.heat(dim=self.dimension, order=self.order)
        elif self.toolboxType == "fluid":
            tb = fluid.fluid(dim=self.dimension)
        else:
            raise ValueError("Unknown toolbox")

        tb.setMesh(mesh)
        tb.setModelProperties(self.model)

        tb.init()

        return tb

    def getField(self, toolbox):
        """Get field of interest from the toolbox

        Args:
            toolbox (Toolbox): Tolbox object

        Raises:
            ValueError: Unknow toolbox

        Returns:
            feelpp_.discr.Element_*: field of the solution
        """
        if self.toolboxType == "heat":
            return toolbox.fieldTemperature()
        elif self.toolboxType == "fluid":
            return toolbox.fieldVelocity()
        else:
            raise ValueError("Unknown toolbox")


    def getToolboxSolution(self, tb, mu):
        """Get the solution of the toolbox tb for the parameter mu

        Args:
            tb (Toolbox): Toolbox object
            mu (parameterSpaceElement): parameter

        Returns:
            feelpp_.discr.Element_*: field of the solution
        """
        assembleToolbox(tb, mu)
        tb.solve()
        return self.getField(tb)


    def createInterpolator(self, domain_tb, image_tb):
        """Create an interpolator between two toolboxes

        Args:
            domain_tb (Toolbox): coarse toolbox
            image_tb  (Toolbox): fine toolbox

        Returns:
            OperatorInterpolation: interpolator object
        """
        if self.toolboxType == "heat":
            Vh_image = image_tb.spaceTemperature()
            Vh_domain = domain_tb.spaceTemperature()
        elif self.toolboxType == "fluid":
            Vh_image = image_tb.spaceVelocity()
            Vh_domain = domain_tb.spaceVelocity()
        interpolator = fi.interpolator(domain = Vh_domain, image = Vh_image, range = image_tb.rangeMeshElements())
        return interpolator

    def getInformations(self):
        """Get information about the model

        Returns:
            dict: information about the model
        """
        info = {}
        info["Ndofs"] = self.Ndofs
        info["Nh"] = self.getFieldSpace().nDof()
        if self.tbCoarse is not None:
            info["NH"] = self.getFieldSpace(True).nDof()
        info["H"] = self.H
        info["h"] = self.h
        info["order"] = self.order
        info["toolboxType"] = self.toolboxType
        info["dimension"] = self.dimension
        if self.outdir is not None:
            info["outdir"] = self.outdir
        return info


class nirbOffline(ToolboxModel):
    """
    Generate the offline part of nirb method 

    Args:
        ToolboxModel (class): class associated to the toolbox model 
    """

    def __init__(self, method="POD", doRectification=True, initCoarse=False, **kwargs) -> None:
        """Initialize the NIRB class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            finemesh_path (str): path to fine mesh file (if geo file, this will be the same as coarse mesh file)
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
            initCoarse (bool, optional): initialize the coarse toolbox. Defaults to False.
            doBiorthonormal (bool, optional): get bi-orthonormalization of reduced basis. Defaults to False.
        """

        super().__init__(**kwargs)

        assert method in ["POD", "Greedy"]
        self.method = method
        self.doRectification = doRectification
        # self.doBiorthonormal = doBiorthonormal

        self.l2ScalarProductMatrix = None
        self.h1ScalarProductMatrix = None
        self.l2ProductBasis = [] # list containing the vector given the column of (l2ScalarProductMatrix @ reeducedBasis) 
        self.reducedBasis = None # list containing the vector of reduced basis function
        self.N = 0 # number of modes

        if self.doRectification or initCoarse:
            super().initCoarseToolbox()
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Initialization done")



    def BiOrthonormalization(self):
        """Bi-orthonormalization of reduced basis
        """
        K = np.zeros((self.N,self.N))
        M = K.copy()

        for i in range(self.N):
            for j in range(self.N):
                K[i,j] = self.h1ScalarProductMatrix.energy(self.reducedBasis[i],self.reducedBasis[j])
                M[i,j] = self.l2ScalarProductMatrix.energy(self.reducedBasis[i],self.reducedBasis[j])

        from scipy import linalg

        eval,evec=linalg.eigh(a=K, b=M, overwrite_a=True, overwrite_b=True) #eigenvalues
        eigenValues = eval.real
        eigenVectors = evec 
        idx = eigenValues.argsort()[::-1]
        eigenValues = eigenValues[idx]
        eigenVectors = evec[:, idx]

        for i in range(self.N):
            eigenVectors[:,i] /= np.sqrt(np.abs(eigenValues[i]))

        oldbasis = self.reducedBasis.copy()
        self.reducedBasis = []

        for i in range(self.N):
            vec = self.Xh.element()
            vec.setZero()
            for j in range(self.N):
                vec = vec + float(eigenVectors[j,i])*oldbasis[j]
                
            self.reducedBasis.append(vec)


    def initProblem(self, numberOfInitSnapshots, Xi_train=None, samplingMode="log-random", computeCoarse=False):
        """Initialize the problem

        Args:
        -----
            numberOfInitSnapshots (int): number of snapshots to use for the initialization
            Xi_train (list of ParameteœrSpaceElement, optional): Train set for algorithm. If None is given, a set of size Ntrain is generated. Defaults to None.
            samplingMode (str, optional): sampling mode in the parameter space.(random, log-random, log-equidistribute, equidistribute) Defaults to "log-random".
            computeCoarse (bool, optional): compute snapshots for coarse toolbox, used for rectification. Defaults to False.
        
        Returns:
        --------
            Xi_train (list of ParameterSpaceElement): parameters used for do build the basis
        """
        if self.doRectification:
            computeCoarse=True

        self.fineSnapShotList = []
        self.coarseSnapShotList = []

        if Xi_train is None:
            s = self.Dmu.sampling()
            s.sampling(numberOfInitSnapshots, samplingMode)
            vector_mu = s.getVector()
        else:
            vector_mu = choices(Xi_train, k=numberOfInitSnapshots)

        if computeCoarse:
            assert self.tbCoarse is not None, f"Coarse toolbox needed for computing coarse Snapshot. set doRectification->True"
            for mu in vector_mu:
                if self.worldcomm.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(self.getToolboxSolution(self.tbFine, mu))
                self.coarseSnapShotList.append(self.getToolboxSolution(self.tbCoarse, mu))

        else:
            for mu in vector_mu:
                if self.worldcomm.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(self.getToolboxSolution(self.tbFine, mu))

        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Number of snapshot computed : {len(self.fineSnapShotList)}" )

        return vector_mu

    def getReducedSolution(self, coarseSolutions, mu, N):
        """Computed the reduced solution for a given parameter and a size of basis

        Args:
        -----
            coarseSolutions (dict of Feel++ solutions): list of coarse solution
            mu (parameterSpaceElement): parameter
            N (int): size of the sub-basis

        Returns:
        --------
            PETSc.Vec: reduced solution u_H^N
        """
        coarseSol = coarseSolutions[mu]
        uHN = np.zeros(N)

        for i in range(N):
            # uHN[i] = self.coarseSnapShotList[i].to_petsc().vec().dot( coarseSol.to_petsc().vec() )
            uHN[i] = self.l2ScalarProductMatrixCoarse.energy( self.coarseSnapShotList[i], coarseSol )
        return uHN


    def initProblemGreedy(self, Ntrain, eps, Xi_train=None, Nmax=50, samplingMode="log-random", computeCoarse=False):
        """Initialize the problem, using a greedy loop

        Args:
        -----
            Ntrain (int): size of the train parameter set
            eps (float): precision of the greedy loop
            Xi_train (list of ParameterSpaceElement, optional): Train set for algorithm. If None is given, a set of size Ntrain is generated. Defaults to None.
            Nmax (int, optional): maximal number of iterations. Defaults to 500.
            samplingMode (str, optional): sampling mode. Defaults to "log-random".
            computeCoarse (bool, optional): compute snapshots for coarse toolbox, used for rectification. Defaults to False.

        Return:
        ------
            a tuple res where
            - res[0] contains the sample of paramters selected
            - res[1] contains a copy of Xi_train sample
            - res[2] contains the evolution of Delta_max

        Raises:
            Exception: Coarse toolbox has not been initialized
        """
        if self.tbCoarse is None:
            raise Exception("Coarse toolbox needed for computing coarse Snapshot. set initCoarse->True in initialization")
        Nmax = min(Nmax, Ntrain)
        if Xi_train is None:
            s = self.Dmu.sampling()
            s.sampling(Ntrain, samplingMode)
            Xi_train = s.getVector()
        Xi_train_copy = Xi_train.copy()

        Delta_star = eps+1
        Deltas_conv = []
        S = []
        self.fineSnapShotList = []
        self.coarseSnapShotList = []
        self.coarseSnapShotListGreedy = []
        N = 0

        # Computation of coarse solutions
        coarseSolutions = {}
        for mu in tqdm(Xi_train_copy, desc="[NIRB] Computing coarse solutions", ncols=120):
            coarseSolutions[mu] = self.getToolboxSolution(self.tbCoarse, mu)

        for i in range(2):
            mu0 = self.Dmu.element()
            S.append(mu0)
            self.fineSnapShotList.append(   self.getToolboxSolution(self.tbFine  , mu0) )
            self.coarseSnapShotList.append( self.getToolboxSolution(self.tbCoarse, mu0) )
            self.coarseSnapShotListGreedy.append(self.getToolboxSolution(self.tbCoarse, mu0) )
            N += 1
            self.orthonormalizeMatL2(self.coarseSnapShotListGreedy)
        Delta_0 = np.abs(self.getReducedSolution(coarseSolutions, Xi_train_copy[0], 2).mean()
                         - self.getReducedSolution(coarseSolutions, Xi_train_copy[1], 1).mean())

        while Delta_star > eps and N < Nmax:
            M = N - 1

            Delta_star = -float('inf')

            for i, mu in enumerate(tqdm(Xi_train_copy,desc=f"[NIRB] Greedy selection", ascii=False, ncols=120)):
                uHN = self.getReducedSolution(coarseSolutions, mu, N)
                uHM = self.getReducedSolution(coarseSolutions, mu, M)

                Delta = np.abs( uHN.mean() - uHM.mean() ) / Delta_0

                if Delta > Delta_star:
                    Delta_star = Delta
                    mu_star = mu
                    idx = i

            S.append(mu_star)
            del Xi_train_copy[idx]
            self.fineSnapShotList.append(  self.getToolboxSolution(self.tbFine  , mu_star))
            self.coarseSnapShotList.append(self.getToolboxSolution(self.tbCoarse, mu_star))
            self.coarseSnapShotListGreedy.append(self.getToolboxSolution(self.tbCoarse, mu_star))
            N += 1
            self.orthonormalizeMatL2(self.coarseSnapShotListGreedy)

            Deltas_conv.append(Delta_star)
            if self.worldcomm.isMasterRank():
                print(f"[nirb] Adding snapshot with mu = {mu_star}")
                print(f"[nirb] Greedy loop done. N = {N}, Delta_star = {Delta_star}")

        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Number of snapshot computed : {N}")

        print(Deltas_conv)

        return S, Xi_train, Deltas_conv


    def generateOperators(self, coarse=False):
        """Assemble L2 and H1 operators, stored in self.l2ScalarProduct and self.h1ScalarProduct
        """
        if self.l2ScalarProductMatrix is None or self.h1ScalarProductMatrix is None:
            # Vh = feelpp.functionSpace(mesh=self.tbFine.mesh(), order=self.order)
            self.l2ScalarProductMatrix = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            self.h1ScalarProductMatrix = FppOp.stiffness(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            self.l2ScalarProductMatrix.to_petsc().mat().assemble()
            self.h1ScalarProductMatrix.to_petsc().mat().assemble()

        if coarse:
            self.l2ScalarProductMatrixCoarse = FppOp.mass(test=self.XH, trial=self.XH, range=feelpp.elements(self.tbCoarse.mesh()))
            self.h1ScalarProductMatrixCoarse = FppOp.stiffness(test=self.XH, trial=self.XH, range=feelpp.elements(self.tbCoarse.mesh()))
            self.l2ScalarProductMatrixCoarse.to_petsc().mat().assemble()
            self.h1ScalarProductMatrixCoarse.to_petsc().mat().assemble()

    def generateReducedBasis(self, tolerance=1.e-6, regulParam=1.e-10):
        """Generate the reduced basis, and store it in the list self.reducedBasis

        Args :
            tolerance(float), optional : the tolerance value for
            regulParam(float), optional : the regularization parameter for rectification
        """
        self.reducedBasis = self.PODReducedBasis(tolerance=tolerance)
        self.orthonormalizeL2()
        self.N = len(self.reducedBasis)
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Number of modes : {self.N}")
        
        if len(self.l2ProductBasis)==0:
            self.getl2ProductBasis()
        # if self.doBiorthonormal:
        #     self.BiOrthonormalization()
        if self.doRectification:
            self.coeffCoarse, self.coeffFine = self.coeffRectification()

    def getl2ProductBasis(self):
        """get the L2 scalar product matrix with reduced basis function 

        """
        assert self.reducedBasis is not None, f"reduced Basis have to be computed before"  

        backend = feelpp.backend(worldcomm=feelpp.Environment.worldCommPtr())

        self.l2ProductBasis = []
        for i in range(self.N):
            vec = backend.newVector(dm=self.Xh.mapPtr()) # beacause the modification of vec will modify the referance at each iteration
            vec.setZero()
            vec.addVector(self.reducedBasis[i], self.l2ScalarProductMatrix)
            self.l2ProductBasis.append(vec)        


    def PODReducedBasis(self, tolerance=1.e-6):
        """
        Computes the reducedOrderBasis using the POD algorithm, from a given list of snapshots contained in self.fineSnapShotList

        Parameters
        ----------
            tolerance (float) : tolerance of the eigen value problem
            target accuracy of the data compression

        Returns
        -------
        ReducedBasis (list) : the reduced basis, of size numberOfModes
        """

        Nsnap = len(self.fineSnapShotList)

        correlationMatrix = PETSc.Mat().createDense(size=Nsnap, comm=MPI.COMM_SELF) # create a dense matrix in sequential mode 
        correlationMatrix.setFromOptions()
        correlationMatrix.setUp()

        for i, snap1 in enumerate(self.fineSnapShotList):
            for j, snap2 in enumerate(self.fineSnapShotList):
                    correlationMatrix[i,j] = self.l2ScalarProductMatrix.energy(snap1,snap2)

        correlationMatrix.assemble()
        eigenValues, eigenVectors =  TruncatedEigenV(correlationMatrix, tolerance)

        Nmode = len(eigenVectors)
        
        for i in range(Nmode):
            eigenVectors[i].scale(1./math.sqrt(abs(eigenValues[i])))

        reducedBasis = []

        for i in range(Nmode):
            vec = self.Xh.element()
            vec.setZero()
            for j in range(Nsnap):
                val = eigenVectors[i].getValue(j)
                vec = vec + val*self.fineSnapShotList[j]

            reducedBasis.append(vec)

        return reducedBasis


    def coeffRectification(self):
        """ Compute the two matrixes used to get rectification matrix R given by :
                
                B_h[i,j] = <U_h(s_i),phi_j >
                B_H[i,j] = <U_H(s_i),phi_j >

        Returns :
                B_h and B_H : the matrix of L2 scalar product between basis function and 
                                fine and coarse snapshot 
        """
        assert len(self.reducedBasis) !=0, f"need computation of reduced basis"

        interpolateOperator = self.createInterpolator(self.tbCoarse, self.tbFine)
        InterpCoarseSnaps = []
        for snap in self.coarseSnapShotList:
            InterpCoarseSnaps.append(interpolateOperator.interpolate(snap))

        coeffCoarse = np.zeros((self.N, self.N))
        coeffFine = np.zeros((self.N, self.N))

        for i in range(self.N):
            for j in range(self.N):
                coeffCoarse[i,j] = self.l2ScalarProductMatrix.energy(InterpCoarseSnaps[i],self.reducedBasis[j])
                coeffFine[i,j] = self.l2ScalarProductMatrix.energy(self.fineSnapShotList[i],self.reducedBasis[j])

        return coeffCoarse, coeffFine

    """
    Check convergence 
    """
    def checkConvergence(self,Ns=10, Xi_test=None, regulParam=1.e-10):
        """
            check convergence of the offline step  
        """
        assert self.tbCoarse is not None, f"Coarse toolbox needed for computing coarse Snapshot. set doRectification->True"

        
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Compute offline convergence error")

        if Xi_test is None:
            s = self.Dmu.sampling()
            s.sampling(Ns,'log-random')
            vector_mu = s.getVector()
        else:
            vector_mu = Xi_test

        Ntest = len(vector_mu)

        interpolationOperator = self.createInterpolator(self.tbCoarse, self.tbFine)

        soltest =[]
        soltestFine=[]
        for mu in vector_mu:
            soltestFine.append(self.getToolboxSolution(self.tbFine, mu))
            uH = self.getToolboxSolution(self.tbCoarse, mu)
            soltest.append(interpolationOperator.interpolate(uH))

        tabcoef = np.zeros((Ntest, self.N))
        tabcoefuh = np.zeros((Ntest, self.N))

        for j in range(Ntest):
            for k in range(self.N):
                tabcoef[j,k] = self.l2ScalarProductMatrix.energy(soltest[j],self.reducedBasis[k])
                tabcoefuh[j,k] = self.l2ScalarProductMatrix.energy(soltestFine[j],self.reducedBasis[k])

        
        Unirb = self.Xh.element()
        UnirbRect = self.Xh.element()
        Uhn = self.Xh.element()

        Error = {'N':[], 'l2(uh-uHn)':[], 'l2(uh-uHn)rec':[], 'l2(uh-uhn)' : [], 'l2(uh-uH)':[]}

        nb = 0
        for i in tqdm(range(1,self.N+1,min(5,self.N)), desc=f"[NIRB] Compute convergence error:", ascii=False, ncols=100):
            nb = i
            # Get rectification matrix
            BH = self.coeffCoarse[:nb, :nb]
            Bh = self.coeffFine[:nb, :nb]
            R = np.zeros((nb,nb)) 
            for s in range(nb):
                R[s,:]=(np.linalg.inv(BH.transpose()@BH+regulParam*np.eye(nb))@BH.transpose()@Bh[:,s])

            for j in range(Ntest):
                Unirb.setZero()
                UnirbRect.setZero()
                Uhn.setZero()
                tabRect = R@tabcoef[j,:nb]

                for k in range(nb): # get reduced sol in a basis function space 
                    Unirb = Unirb + float(tabcoef[j,k])*self.reducedBasis[k]
                    UnirbRect = UnirbRect + float(tabRect[k])*self.reducedBasis[k]
                    Uhn = Uhn + float(tabcoefuh[j,k])*self.reducedBasis[k]

                Unirb = Unirb - soltestFine[j]
                UnirbRect = UnirbRect - soltestFine[j]
                Uhint = soltestFine[j] - soltest[j]
                Uhn = Uhn - soltestFine[j]
                err1 = np.sqrt(abs(self.l2ScalarProductMatrix.energy(Unirb,Unirb)))
                err2 = np.sqrt(abs(self.l2ScalarProductMatrix.energy(UnirbRect,UnirbRect))) 
                err3 = np.sqrt(abs(self.l2ScalarProductMatrix.energy(Uhn,Uhn)))  
                err4 = np.sqrt(abs(self.l2ScalarProductMatrix.energy(Uhint,Uhint))) 

                Error['l2(uh-uHn)'].append(err1)
                Error['l2(uh-uHn)rec'].append(err2)
                Error['l2(uh-uhn)'].append(err3)
                Error['l2(uh-uH)'].append(err4)
                Error["N"].append(nb)

        return Error
    
    """
    Handle Gram-Schmidt orthogonalization
    """
    def orthonormalizeH1(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using H1 norm
        (the optional argument is not needed)
        """
        Z = self.reducedBasis

        Z[0] = Z[0] * (1./math.sqrt(abs(self.h1ScalarProductMatrix.energy(Z[0],Z[0]))))

        for n in range(1, self.N):
            s = self.Xh.element()
            s.setZero()
            for m in range(n):
                s = s + self.h1ScalarProductMatrix.energy(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp * (1./math.sqrt(abs(self.h1ScalarProductMatrix.energy(z_tmp,z_tmp))))

        self.reducedBasis = Z

        if not (self.checkH1Orthonormalized() ) and nb < 10:
            self.orthonormalizeH1(nb=nb+1)
        elif self.worldcomm.isMasterRank():
            print(f"[NIRB] Gram-Schmidt H1 orthonormalization done after {nb+1} step"+['','s'][nb>0])


    def orthonormalizeL2(self, nb=0, tol=1.e-6):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using L2 norm
        (the optional argument is not needed)

        Args : 
            nb (int) : number of orthonormalization to run
            tol (float) : tolerence 
        """

        Z = self.reducedBasis
        Z[0] = Z[0] * (1./math.sqrt(abs(self.l2ScalarProductMatrix.energy(Z[0],Z[0]))))

        for n in range(1, self.N):
            s = self.Xh.element()
            s.setZero()
            for m in range(n):
                s = s + self.l2ScalarProductMatrix.energy(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp * (1./math.sqrt(abs(self.l2ScalarProductMatrix.energy(z_tmp,z_tmp))))

        self.reducedBasis = Z

        if not (self.checkL2Orthonormalized(tol=tol)) and nb < 10:
            self.orthonormalizeL2(nb=nb+1)
        elif self.worldcomm.isMasterRank():
            print(f"[NIRB] Gram-Schmidt L2 orthonormalization done after {nb+1} step"+['','s'][nb>0])

    def orthonormalizeMatL2(self, Z):
        """Ortohnormalize the matrix Z using L2 norm

        Args:
            Z (list of feelpp._discr.Element): list of feelpp functions

        Returns:
            Z: the matrix orthonormalized, inplace
        """
        N = len(Z)

        if N == 1:
            normS = self.l2ScalarProductMatrixCoarse.energy(Z[-1], Z[-1])
            Z[-1] = 1 / np.sqrt(normS) * Z[-1]
        else:
            s = self.XH.element()
            s.setZero()
            for m in range(N):
                Z[-1] =  Z[-1] - self.l2ScalarProductMatrixCoarse.energy(Z[-1], Z[m]) * Z[m]
            normS = self.l2ScalarProductMatrixCoarse.energy(Z[-1], Z[-1])
            Z[-1] = 1 / np.sqrt(normS) * Z[-1]
        return Z

    def checkH1Orthonormalized(self, tol=1e-8):
        """Check if the reduced basis is H1 orthonormalized.

        Args:
            tol (float, optional): Tolerance. Defaults to 1e-8.

        Returns:
            bool: True if the reduced basis is H1 orthonormalized
        """
        
        matH1 = np.zeros((self.N,self.N))

        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    if abs(matH1[i, j] - 1) > tol:
                        print(f"[NIRB] not H1 ortho : pos = {i} {j} val = {abs(matH1[i,j]-1.)} tol = {tol}")
                        return False
                else:
                    if abs(matH1[i, j]) > tol :
                        print(f"[NIRB] not H1 ortho : pos = {i} {j} val = {abs(matH1[i,j])} tol = {tol}")
                        return False
        return True

    def checkL2Orthonormalized(self, tol=1e-8):
        """Check if the reduced basis is L2 orthonormalized.

        Args:
            tol (float, optional): Tolerance. Defaults to 1e-8.

        Returns:
            bool: True if the reduced basis is L2 orthonormalized
        """ 

        matL2 = np.zeros((self.N,self.N))

        for i in range(self.N):
            for j in range(self.N):
                matL2[i,j] = self.l2ScalarProductMatrix.energy(self.reducedBasis[i], self.reducedBasis[j])
                if i == j:
                    if abs(matL2[i, j] - 1) > tol:
                        print(f"[NIRB] not L2 ortho : pos = [{i} , {j}] val = {abs(matL2[i,j] - 1.)} tol = {tol}")
                        return False
                else:
                    if abs(matL2[i, j]) > tol :
                        print(f"[NIRB] not L2 ortho : pos = [{i} , {j}] val = {abs(matL2[i,j])} tol = {tol}")
                        return False
        return True

    def saveData(self, path="./", force=False):
        """ 
        Save the data generated by the offline phase

        Args :
            path (str, optional): Path where files are saved. Defaults to "./".
            force (bool, optional): Force saving, even if files are already present. Defaults to False.
        """

        reducedPath = os.path.join(path, 'reducedBasis')
        reducedFilename = 'reducedBasis'

        l2productPath = os.path.join(path,  'l2productBasis')
        l2productFilename = 'l2productBasis'

        if self.worldcomm.isMasterRank():
            if os.path.isdir(path) and not force:
                print(f"[NIRB] Directory {path} already exists. Rerun with force=True to force saving")
                return
            if not os.path.isdir(path):
                os.makedirs(path)
            if not os.path.isdir(reducedPath):
                os.makedirs(reducedPath)
            if not os.path.isdir(l2productPath):
                os.makedirs(l2productPath)

        self.worldcomm.globalComm().Barrier()

        for i in range(len(self.reducedBasis)):
            self.reducedBasis[i].save(reducedPath, reducedFilename, suffix=str(i))

        for i in range(len(self.l2ProductBasis)):
            vec = self.Xh.element(self.l2ProductBasis[i])
            vec.save(l2productPath, l2productFilename, suffix=str(i))

        nploc = self.worldcomm.localRank()
        coeffCoarseFile = os.path.join(path,f"coeffcoarse_np{nploc}")
        coeffFineFile = os.path.join(path,f"coefffine_np{nploc}")
        if self.doRectification:
            # if self.worldcomm.isMasterRank():
            np.save(coeffCoarseFile, self.coeffCoarse)
            np.save(coeffFineFile, self.coeffFine)

        self.outdir = os.path.abspath(path)
        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Data saved in {self.outdir}")

### ONLINE PHASE ###

class nirbOnline(ToolboxModel):
    """a class to generate the online part of nirb method 

    Args:
        ToolboxModel (class): class associated to the toolbox model 
    """

    def __init__(self, **kwargs):
        """Initialize the NIRB online class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            finemesh_path (str): path to fine mesh file (if geo file, this will be the same as coarse mesh file)
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """

        super().__init__(**kwargs)

        self.doRectification = kwargs['doRectification']

        self.l2ProductBasis = None
        self.reducedBasis = None
        self.N = 0

        super().initCoarseToolbox()

        self.interpolationOperator = self.createInterpolator(self.tbCoarse, self.tbFine)
        self.exporter = None

        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Initialization done")

    def getCompressedSol(self, mu=None, solution=None):
        """
        get the projection of given solution from fine mesh into the reduced space

        Parameters
        ----------
        solution (feelpp._discr.Element) : the solution to be projected
        mu (ParameterSpaceElement) : parameter to compute the solution if not given

        return :
        --------
        compressedSol (numpy.ndarray) : the compressed solution, of size (self.N)
        """
        assert (mu != None) or (solution != None), f"One of the arguments must be given: solution or mu"

        if solution is None:
            sol = self.getInterpSol(mu)
        else :
            sol = solution

        compressedSol = np.zeros(self.N)

        for i in range(self.N):
            compressedSol[i] = self.l2ProductBasis[i].to_petsc().dot(sol.to_petsc())

        return compressedSol

    def getInterpSol(self, mu):
        """Get the interpolated solution from coarse mesh to fine one

        Args:
        -----
            mu (ParameterSpaceElement): parameter

        Returns:
        --------
            interpSol (feelpp._discr.Element): interpolated solution on fine mesh
        """
        interpSol = self.solveOnline(mu)[1]
        return interpSol

    def getOnlineSol(self,mu):
        """Get the Online nirb approximate solution
        """

        onlineSol = self.Xh.element()
        onlineSol.setZero()

        compressedSol = self.getCompressedSol(mu)

        if self.doRectification:
            coef = self.RectificationMat@compressedSol
            for i in range(self.N):
                onlineSol = onlineSol + float(coef[i])*self.reducedBasis[i]
            
            if self.worldcomm.isMasterRank():
                print("[NIRB] Solution computed with Rectification post-process ")
        else :
            for i in range(self.N):
                onlineSol = onlineSol + float(compressedSol[i])*self.reducedBasis[i]

        # to export, call self.exportField(onlineSol, "U_nirb")


        return onlineSol

    def getRectification(self, coeffCoarse, coeffFine, lambd=1.e-10):
        """compute rectification matrix associated to number of basis function selected 
                in the online step.

        Args:
            coeffCoarse (numpy.array): the coefficient matrix gieven by (U^H_i, \phi_j)
            coeffFine (numpy.array): the coefficient matrix (U^h_i,\phi_j)
            lambd (float, optional): regularization parameter. Defaults to 1.e-10.
        """
        BH = coeffCoarse[:self.N, :self.N]
        Bh = coeffFine[:self.N, :self.N]
        R = np.zeros((self.N,self.N))

        #Thikonov regularization (AT@A +lambda I_d)^-1@(AT@B)
        for i in range(self.N):
            R[i,:]=(np.linalg.inv(BH.transpose()@BH+lambd*np.eye(self.N))@BH.transpose()@Bh[:,i])

        return R

    def generateOperators(self,h1=False):
        """Assemble L2 and H1 operators associated to the fine toolbox 
        """
        if h1:
            l2Mat = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            h1Mat = FppOp.stiffness(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            l2Mat.to_petsc().mat().assemble()
            h1Mat.to_petsc().mat().assemble()
            return l2Mat, h1Mat
        else :
            l2Mat = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            l2Mat.to_petsc().mat().assemble()
            return l2Mat

    def loadData(self, nbSnap=None, path="./", regulParam=1.e-10):
        """ 
        Load the data generated by the offline phase

        Args : 
            nbSnap (int, optional) : number of basis function. If not given, all the basis are loaded
            path (str, optional): Path where files are saved. Defaults to "./".
            regulParam (float, optional): regularization parameter of the rectification. Defaults to 1.e-10
        Returns :
            int: error code, 0 if all went well, 1 if not
        """
        print(f"[NIRB] Loading data from {os.path.abspath(path)}")

        reducedPath = os.path.join(path,'reducedBasis')
        reducedFilename = 'reducedBasis'

        l2productPath = os.path.join(path, 'l2productBasis')
        l2productFilename = 'l2productBasis'
        
        if self.worldcomm.isMasterRank():
            if not os.path.isdir(path):
                print(f"[NIRB] Error : Could not find path {path}")
                return 1
            if not os.path.isdir(reducedPath):
                print(f"[NIRB] Error : Could not find path {reducedPath}")
                return 1
            if not os.path.isdir(l2productPath):
                print(f"[NIRB] Error : Could not find path {l2productPath}")
                return 1

        self.worldcomm.globalComm().Barrier()

        self.reducedBasis = []
        self.l2ProductBasis = []

        if nbSnap is None:
            import glob
            Nreduce = len(glob.glob(os.path.join(reducedPath, "*.h5")))
            Nl2 = len(glob.glob(os.path.join(l2productPath, "*.h5")))

            assert Nreduce == Nl2, f"different number of files, {Nreduce} != {Nl2}"
            nbSnap = Nreduce

        for i in range(nbSnap):
            vec = self.Xh.element()
            vec.load(reducedPath, reducedFilename, suffix=str(i))
            self.reducedBasis.append(vec)
        
        self.N = nbSnap

        for i in range(nbSnap):
            vec = self.Xh.element()
            vec.load(l2productPath, l2productFilename, suffix=str(i))
            self.l2ProductBasis.append(vec)

        # coeffCoarseFile = os.path.join(path,'coeffcoarse.npy')
        # coeffFineFile = os.path.join(path,'coefffine.npy')
        nploc = self.worldcomm.localRank()
        coeffCoarseFile = os.path.join(path,f"coeffcoarse_np{nploc}.npy")
        coeffFineFile = os.path.join(path,f"coefffine_np{nploc}.npy")
        if self.doRectification:
            coeffCoarse = np.load(coeffCoarseFile)
            coeffFine = np.load(coeffFineFile)
            self.RectificationMat = self.getRectification(coeffCoarse, coeffFine, lambd=regulParam)

        if self.worldcomm.isMasterRank():
            print(f"[NIRB] Data loaded from {os.path.abspath(path)}")
            print(f"[NIRB] Number of basis functions loaded : {self.N}")
        
        return 0

    def normMat(self, Mat, u):
        """ Compute the norm associated to matrix Mat

        Args:
            Mat (feelpp.__alg.SparseMatrix): the matrix 
            u (feelpp.__discr_Vector): the vector to compute the norm  
        """
        return np.sqrt(np.abs(Mat.energy(u,u)))

    def solveOnline(self, mu):
        """Retrun the interpolated FE-solution from coarse mesh to fine one u_{Hh}^calN(mu)
            Solve in Coarse mesh and interpolate in fine mesh

        Args:
            mu (ParameterSpaceElement): parameter.

        Returns:
            feelpp._discr.Element: FE solution on coarse mesh
            feelpp._discr.Element: interpolated solution on fine mesh
        """
        if self.tbCoarse is None:
            super().initCoarseToolbox()

        coarseSol = self.getToolboxSolution(self.tbCoarse, mu)
        interpolatedSol = self.interpolationOperator.interpolate(coarseSol)

        # Export
        # call self.exportField(interpolatedSol, "U_interp")

        return coarseSol, interpolatedSol


    def initExporter(self, name, toolbox="fine"):
        """init feelpp exporter

        Args:
            name (str): name of the exporter
            toolbox (str, optional): mesh to use, "fine" or "coarse". Defaults to "fine".
        """
        if toolbox == "fine":
            self.exporter = feelpp.exporter(self.tbFine.mesh(), name)
        elif toolbox == "coarse":
            self.exporter = feelpp.exporter(self.tbCoarse.mesh(), name)

    def exportField(self, field, name):
        """export a field to the exporter, if it has already been initialized

        Args:
            field (feelpp._discr.Element): field to export
            name (str): name of the field
        """
        if self.exporter is not None:
            if self.order == 1:
                self.exporter.addP1c(name, field)
            elif self.order == 2:
                self.exporter.addP2c(name, field)
        else:
            print("Exporter not initialized, please call initExporter() first")

    def saveExporter(self):
        """save the exporter to disk
        """
        if self.exporter is not None:
            self.exporter.save()
        else:
            print("Exporter not initialized, please call initExporter() first")
