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
import math 


import os

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class ToolboxModel():

    def __init__(self, dim, H, h, toolboxType, model_path, geo_path, order=1, **kwargs) -> None:
        """Initialize the toolbox model class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            geo_path (str): path to geo file
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
        self.geo_path = feelpp.Environment.expand(geo_path)

        self.tbCoarse  = None
        self.tbFine    = None

        self.Xh = None

        self.initModel()
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Initialization done")


    def initModel(self):
        """Initialize the model
        """
        self.model = feelpp.readJson(self.model_path)
        self.tbFine = self.setToolbox(self.h)
        self.Xh = feelpp.functionSpace(mesh=self.tbFine.mesh(), order=self.order)
        self.Dmu = loadParameterSpace(self.model_path)
        self.Ndofs = self.getFieldSpace().nDof()

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of nodes on the fine mesh : {self.Ndofs}")

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
        self.tbCoarse = self.setToolbox(self.H)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of nodes on the coarse mesh : {self.tbCoarse.mesh().numGlobalPoints()}")


    def setToolbox(self, hsize):
        """Set up the toolbox object for the given model and mesh

        Args:
            hsize (float): size of the mesh

        Returns:
            Toolbox: toolbox object
        """

        # load meshes
        mesh_ = feelpp.mesh(dim=self.dimension, realdim=self.dimension)
        mesh = feelpp.load(mesh_, self.geo_path, hsize)

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
        return info


class nirbOffline(ToolboxModel):

    def __init__(self, method="POD", doRectification=True, **kwargs) -> None:
        """Initialize the NIRB class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            geo_path (str): path to geo file
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
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

        if self.doRectification:
            super().initCoarseToolbox()
        if feelpp.Environment.isMasterRank():
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



    def initProblem(self, numberOfInitSnapshots, samplingMode="log-random", computeCoarse=False):
        """Initialize the problem

        Args:
            numberOfInitSnapshots (int): number of snapshots to use for the initialization
            samplingMode (str, optional): sampling mode in the parameter space.(random, log-random, log-equidistribute, equidistribute) Defaults to "log-random".
            computeCoarse (bool, optional): compute snapshots for coarse toolbox, used for rectification. Defaults to False.
        """
        if self.doRectification:
            computeCoarse=True

        self.fineSnapShotList = []
        self.coarseSnapShotList = []
 
        s = self.Dmu.sampling()
        s.sampling(numberOfInitSnapshots, samplingMode)
        vector_mu = s.getVector()
        # vector_mu = samplingEqui(self.model_path, numberOfInitSnapshots)
        # vector_mu = comm.Bcast(vector_mu, root=0)

        if computeCoarse:
            assert self.tbCoarse is not None, f"Coarse toolbox needed for computing coarse Snapshot. set doRectification->True"
            for mu in vector_mu:
                if feelpp.Environment.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(self.getToolboxSolution(self.tbFine, mu))
                self.coarseSnapShotList.append(self.getToolboxSolution(self.tbCoarse, mu))

        else:
            for mu in vector_mu:
                if feelpp.Environment.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(self.getToolboxSolution(self.tbFine, mu))

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of snapshot computed : {len(self.fineSnapShotList)}" )


    def generateOperators(self):
        """Assemble L2 and H1 operators, stored in self.l2ScalarProduct and self.h1ScalarProduct
        """
        if self.l2ScalarProductMatrix is None or self.h1ScalarProductMatrix is None:
            # Vh = feelpp.functionSpace(mesh=self.tbFine.mesh(), order=self.order)
            self.l2ScalarProductMatrix = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            self.h1ScalarProductMatrix = FppOp.stiffness(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            self.l2ScalarProductMatrix.to_petsc().mat().assemble()
            self.h1ScalarProductMatrix.to_petsc().mat().assemble()

    def generateReducedBasis(self, tolerance=1.e-6, regulParam=1.e-10):
        """Generate the reduced basis, and store it in the list self.reducedBasis

        Args :
            tolerance(float), optional : the tolerance value for
            regulParam(float), optional : the regularization parameter for rectification
        """
        self.reducedBasis = self.PODReducedBasis(tolerance=tolerance)
        self.N = len(self.reducedBasis)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of modes : {self.N}")
        
        if len(self.l2ProductBasis)==0:
            self.getl2ProductBasis()
        # if self.doBiorthonormal:
        #     self.BiOrthonormalization()
        if self.doRectification:
            self.RectificationMat = self.Rectification(lambd=regulParam)

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

    def Rectification(self, lambd=1e-10):
        """ Compute the rectification matrix R given by :
                R = B_h*(B_H)^-1
                with B_h[i,j] = <U_h(s_i),\phi_j >
                and B_H[i,j] = <U_H(s_i),\phi_j >

        Args :
            lambd (float) : Tikonov regularization parameter

        Returns :
            R (petsc.Mat) : the rectification matrix
        """
        assert len(self.reducedBasis) !=0, f"need computation of reduced basis"

        interpolateOperator = self.createInterpolator(self.tbCoarse, self.tbFine)
        InterpCoarseSnaps = []
        for snap in self.coarseSnapShotList:
            InterpCoarseSnaps.append(interpolateOperator.interpolate(snap))

        R = np.zeros((self.N,self.N))
        BH = R.copy()
        Bh = R.copy()

        for i in range(self.N):
            for j in range(self.N):
                BH[i,j] = self.l2ScalarProductMatrix.energy(InterpCoarseSnaps[i],self.reducedBasis[j])
                Bh[i,j] = self.l2ScalarProductMatrix.energy(self.fineSnapShotList[i],self.reducedBasis[j])


        #Thikonov regularization (AT@A +lambda I_d)^-1
        for i in range(self.N):
            R[i,:]=(np.linalg.inv(BH.transpose()@BH+lambd*np.eye(self.N))@BH.transpose()@Bh[:,i])

        return R


    """
    Handle Gram-Schmidt orthogonalization
    """
    def orthonormalizeH1(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using H1 norm
        (the optional argument is not needed)
        """
        Z = self.reducedBasis

        Z[0] = Z[0] * (1./self.h1ScalarProductMatrix.energy(Z[0],Z[0]))

        for n in range(1, self.N):
            s = self.Xh.element()
            s.setZero()
            for m in range(n):
                s = s + self.h1ScalarProductMatrix.energy(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp * (1./self.h1ScalarProductMatrix.energy(z_tmp,z_tmp))

        self.reducedBasis = Z

        if not (self.checkH1Orthonormalized() ) and nb < 10:
            self.orthonormalizeH1(nb=nb+1)
        elif rank == 0:
            # pass
            print(f"[NIRB] Gram-Schmidt H1 orthonormalization done after {nb+1} step"+['','s'][nb>0])

    def orthonormalizeL2(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using L2 norm
        (the optional argument is not needed)
        """

        Z = self.reducedBasis

        Z[0] = Z[0] * (1./self.l2ScalarProductMatrix.energy(Z[0],Z[0]))

        for n in range(1, self.N):
            s = self.Xh.element()
            s.setZero()
            for m in range(n):
                s = s + self.l2ScalarProductMatrix.energy(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp * (1./self.l2ScalarProductMatrix.energy(z_tmp,z_tmp))

        self.reducedBasis = Z

        if not (self.checkL2Orthonormalized() ) and nb < 10:
            self.orthonormalizeL2(nb=nb+1)
        elif rank == 0:
            # pass
            print(f"[NIRB] Gram-Schmidt L2 orthonormalization done after {nb+1} step"+['','s'][nb>0])

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
                        return False
                    # assert abs(matH1[i, j] - 1) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"
                else:
                    if abs(matH1[i, j]) > tol :
                        return False
                    # assert abs(matH1[i, j]) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"
        return True

    def checkL2Orthonormalized(self, tol=1e-6):
        """Check if the reduced basis is L2 orthonormalized.

        Args:
            tol (float, optional): Tolerance. Defaults to 1e-8.

        Returns:
            bool: True if the reduced basis is L2 orthonormalized
        """
        
        assert len(self.l2ProductBasis) == len(self.reducedBasis) != 0 

        matL2 = np.zeros((self.N,self.N))

        for i in range(self.N):
            for j in range(self.N):
                matL2[i,j] = self.l2ProductBasis[i].to_petsc().dot(self.reducedBasis[j].to_petsc())
                if i == j:
                    if abs(matL2[i, j] - 1) > tol:
                        return False
                else:
                    if abs(matL2[i, j]) > tol :
                        return False
        return True

    def saveData(self, path="./", force=False):
        """ 
        Save the data generated by the offline phase

        Args : 
            path (str, optional): Path where files are saved. Defaults to "./".
            force (bool, optional): Force saving, even if files are already present. Defaults to False.
        """


        reducedPath = path +'/reducedBasis'
        reducedFilename = 'reducedBasis'

        l2productPath = path + '/l2productBasis'
        l2productFilename = 'l2productBasis'

        if feelpp.Environment.isMasterRank():
            if os.path.isdir(path) and not force:
                print(f"[NIRB] Directory {path} already exists. Rerun with force=True to force saving")
                return
            elif not os.path.isdir(path):
                os.mkdir(path)
            if not os.path.isdir(reducedPath):
                os.mkdir(reducedPath)
            if not os.path.isdir(l2productPath):
                os.mkdir(l2productPath)
        
        comm.Barrier()

        
        for i in range(len(self.reducedBasis)):
            self.reducedBasis[i].save(reducedPath, reducedFilename, suffix=str(i))

        for i in range(len(self.l2ProductBasis)):
            vec = self.Xh.element(self.l2ProductBasis[i])
            vec.save(l2productPath, l2productFilename, suffix=str(i))

        rectificationFile = 'rectification'
        if self.doRectification:
            np.save(rectificationFile, self.RectificationMat)

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Data saved in {os.path.abspath(path)}")

### ONLINE PHASE ###

class nirbOnline(ToolboxModel):

    def __init__(self, **kwargs):
        """Initialize the NIRB online class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            model_path (str): path to json file
            geo_path (str): path to geo file
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

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Initialization done")

    def getCompressedSol(self, mu=None, solution=None):
        """
        get the projection of given solution from fine mesh into the reduced space

        Parameters
        ----------
        solution (feelpp._discr.Element) : the solution to be projected
        mu (ParameterSpaceElement) : parameter to compute the solution if not given 

        return :
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
        """Get the interpolate solution from coarse mesh to fine one

        Args:
            mu (ParameterSpaceElement): parameter

        Returns:
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
            
            if feelpp.Environment.isMasterRank():
                print("[NIRB] Solution computed with Rectification post-process ")
        else :
            for i in range(self.N):
                onlineSol = onlineSol + float(compressedSol[i])*self.reducedBasis[i]

        # to export, call self.exportField(onlineSol, "U_nirb")


        return onlineSol

    def generateOperators(self,h1=False):
        """Assemble L2 and H1 operators associated to the fine toolbox 
        """
        if h1:
            l2Mat = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            h1Mat = FppOp.stiffness(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            l2Mat.mat().assemble()
            h1Mat.mat().assemble()
            return l2Mat, h1Mat
        else :
            l2Mat = FppOp.mass(test=self.Xh, trial=self.Xh, range=feelpp.elements(self.tbFine.mesh()))
            l2Mat.mat().assemble()
            return l2Mat

    def loadData(self, path="./"):
        """ 
        Load the data generated by the offline phase

        Args : 
            path (str, optional): Path where files are saved. Defaults to "./".
        """


        reducedPath = path +'/reducedBasis'
        reducedFilename = 'reducedBasis'

        l2productPath = path + '/l2productBasis'
        l2productFilename = 'l2productBasis'
        
        if feelpp.Environment.isMasterRank():
            if not os.path.isdir(path):
                print(f"[NIRB] Error : Could not find path {path}.")
                return None 
            if not os.path.isdir(reducedPath):
                print(f"[NIRB] Error : Could not find path {reducedPath}.")
                return None 
            if not os.path.isdir(l2productPath):
                print(f"[NIRB] Error : Could not find path {l2productPath}.")
                return None 
                    

        comm.Barrier()

        self.reducedBasis = []
        self.l2ProductBasis = []


        import glob
        Nreduce = len(glob.glob(os.path.join(reducedPath, "*.h5")))
        Nl2 = len(glob.glob(os.path.join(l2productPath, "*.h5")))

        assert Nreduce == Nl2, f"different number of files, {Nreduce} != {Nl2}"

        for i in range(Nreduce):
            vec = self.Xh.element()
            vec.load(reducedPath, reducedFilename, suffix=str(i))
            self.reducedBasis.append(vec)
        
        self.N = len(self.reducedBasis)

        for i in range(Nl2):
            vec = self.Xh.element()
            vec.load(l2productPath, l2productFilename, suffix=str(i))
            self.l2ProductBasis.append(vec)

        assert self.N == len(self.l2ProductBasis), f"{self.N} != {len(self.l2ProductBasis)} "

        rectificationFile = 'rectification.npy'
        if self.doRectification:
            self.RectificationMat = np.load(rectificationFile)

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Data loaded from {os.path.abspath(path)}")
            print(f"[NIRB] Number of basis functions loaded : {self.N}")

    def normMat(self, Mat, u):
        """ Compute the norm associated to matrix Mat

        Args:
            Mat (feelpp.__alg.SparseMatrix): the matrix 
            u (feelpp.__discr_Vector): the vector to compute the norm  
        """
        return np.sqrt(np.abs(Mat.energy(u,u)))

    def solveOnline(self, mu):
        """Retrun the interpolated FE-solution from coarse mesh to fine one u_{Hh}^calN(\mu)
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
            self.exporter = feelpp.Exporter(self.tbFine.mesh(), name)
        elif toolbox == "coarse":
            self.exporter = feelpp.Exporter(self.tbCoarse.mesh(), name)

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