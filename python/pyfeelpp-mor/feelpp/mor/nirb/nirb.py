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


import os

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class ToolboxModel():

    def __init__(self, dimension, H, h, toolboxType, cfg_path, model_path, geo_path, order=1) -> None:
        """Initialize the toolbox model class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            cfg_path (str): path to cfg file
            model_path (str): path to json file
            geo_path (str): path to geo file
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """

        assert dimension in [2,3]
        self.dimension = dimension
        # assert method in ["POD", "Greedy"]
        self.H = H
        self.h = h if h is not None else H**2
        # self.method = method
        self.order = order

        # self.doRectification = doRectification

        self.toolboxType = toolboxType
        assert toolboxType in ["heat", "fluid"], "toolboxType must be 'heat' or 'fluid'"
        self.cfg_path = cfg_path
        self.model_path = model_path
        self.geo_path = geo_path

        self.tbCoarse  = None
        self.tbFine    = None

        self.Xh = None

        self.initModel()
        # if self.doRectification:
            # self.initRectification()
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Initialization done")


    def initModel(self):
        """Initialize the model
        """
        self.model = feelpp.readJson(self.model_path)
        self.tbFine = self.setToolbox(self.h)
        self.Xh = feelpp.functionSpace(mesh=self.tbFine.mesh(), order=self.order)
        self.Dmu = loadParameterSpace(self.model_path)
        self.Ndofs = self.tbFine.mesh().numGlobalPoints()

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of nodes on the fine mesh : {self.Ndofs}")


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


class nirbOffline(ToolboxModel):

    def __init__(self, dimension, H, h, toolboxType, cfg_path, model_path, geo_path,
        method="POD", order=1, doRectification=True) -> None:
        """Initialize the NIRB class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            cfg_path (str): path to cfg file
            model_path (str): path to json file
            geo_path (str): path to geo file
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """

        super().__init__(dimension, H, h, toolboxType, cfg_path, model_path, geo_path, order)

        assert method in ["POD", "Greedy"]
        self.method = method
        self.doRectification = doRectification

        self.l2ScalarProductMatrix = None
        self.h1ScalarProductMatrix = None
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

        # K = PETSc.Mat().createDense(size=(self.N,self.N))
        # K.setFromOptions()
        # K.setUp()
        # K.assemble()
        # M = K.copy()

        Udof, Umode = self.reducedBasis.createVecs()
        xr = Udof.copy()

        for i in range(self.N):
            xr[:] = self.reducedBasis[:,i]
            self.h1ScalarProductMatrix.to_petsc().mat().mult(xr,Udof)
            self.reducedBasis.mult(Udof,Umode)  # ??
            K[i,:] = Umode[:]

            self.l2ScalarProductMatrix.to_petsc().mat().mult(xr,Udof)
            self.reducedBasis.mult(Udof,Umode)
            M[i,:] = Umode[:]

        from scipy import linalg

        vc,vr=linalg.eigh(K, b=M) #eigenvalues
        eigenValues = vc.real
        idx = eigenValues.argsort()[::-1]
        eigenValues = eigenValues[idx]
        eigenVectors = vr[:, idx]

        for i in range(self.N):
            eigenVectors[i,:] /= np.sqrt(eigenValues[i])


        eigenMatrix = PETSc.Mat().createDense([self.N,self.N])
        eigenMatrix.setFromOptions()
        eigenMatrix.setUp()

        eigenMatrix[:,:] = eigenVectors[:,:]

        eigenMatrix.assemble()

        bb = self.reducedBasis.copy()
        eigenMatrix.matMult(bb,self.reducedBasis)

        bb.destroy()
        eigenMatrix.destroy()



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
        """Generate the reduced basis, and store it in self.reducedBasis

        Args :
            tolerance(float), optional : the tolerance value for
            regulParam(float), optional : the regularization parameter for rectification
        """
        self.reducedBasis = self.PODReducedBasis(tolerance=tolerance)
        self.N = self.reducedBasis.size[1]
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of modes : {self.N}")
        if self.doRectification:
            self.RectificationMat = self.Rectification(lambd=regulParam)


    def PODReducedBasis(self, tolerance=1.e-6):
        """
        Computes the reducedOrderBasis using the POD algorithm, from a given list of snapshots contained in self.fineSnapShotList

        Parameters
        ----------
            tolerance (float) : tolerance of the eigen value problem
            target accuracy of the data compression

        Returns
        -------
        ReducedBasis (petsc.Mat) : the reduced basis, of size (numberOfModes, numberOfDOFs)
        """

        Nsnap = len(self.fineSnapShotList)
        correlationMatrix = PETSc.Mat().create()
        correlationMatrix.setSizes([Nsnap,Nsnap])
        correlationMatrix.setFromOptions()
        correlationMatrix.setUp()

        for i, snap1 in enumerate(self.fineSnapShotList):
            for j, snap2 in enumerate(self.fineSnapShotList):
                    correlationMatrix[i,j] = self.scalarL2(snap1.to_petsc().vec(),snap2.to_petsc().vec())

        correlationMatrix.assemble()
        eigenValues, eigenVectors =  TruncatedEigenV(correlationMatrix, tolerance) # truncate only eigenvalue >0

        Nmode = len(eigenVectors)
        for i in range(Nmode):
            eigenVectors[i] /= np.sqrt(np.abs(eigenValues[i]))

        LS = []
        for i in range(Nsnap):
            LS.append(self.fineSnapShotList[i].to_petsc().vec()[:])


        reducedOrderBasis = PETSc.Mat().createDense(size=(self.Ndofs, Nmode), comm=comm)
        reducedOrderBasis.setFromOptions()
        reducedOrderBasis.setUp()
        reducedOrderBasis.assemble()

        # reducedOrderBasis.transpose()

        # vec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        # vec.setSizes(self.Ndofs)
        # vec.setFromOptions()
        # vec.setUp()

        vec = np.zeros(self.Ndofs)

        for i in range(Nmode):
            vec.fill(0.)
            for j in range(Nsnap):
                vec += float(eigenVectors[i][j])*LS[j]

            reducedOrderBasis[:,i] = vec 

            # r = reducedOrderBasis.getDenseColumnVec(i)
            # vec.copy(r)
            # reducedOrderBasis.restoreDenseColumnVec(i)


        reducedOrderBasis.assemble()
        # reducedOrderBasis.transpose()

        return reducedOrderBasis

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
        assert self.N == self.reducedBasis.size[1], f"need computation of reduced basis"

        interpolateOperator = self.createInterpolator(self.tbCoarse, self.tbFine)
        CoarseSnaps = []
        for snap in self.coarseSnapShotList:
            CoarseSnaps.append(interpolateOperator.interpolate(snap))

        BH = np.zeros((self.N,self.N))
        Bh = BH.copy()


        R = PETSc.Mat().createDense(size=(self.N,self.N))
        R.setFromOptions()
        R.setUp()
        R.assemble()

        lfine = []
        lcoarse = []

        Usnap, Udof = self.reducedBasis.createVecs()

        # vec = PETSc.Vec().create(comm=comm)
        # vec.setSizes(self.Ndofs)
        # vec.setFromOptions()
        # vec.setUp()

        # Udof = PETSc.Vec().create(comm=comm)
        # Udof.setSizes(self.Ndofs)
        # Udof.setFromOptions()
        # Udof.setUp()

        # Usnap = PETSc.Vec().create(comm=comm)
        # Usnap.setSizes(self.N)
        # Usnap.setFromOptions()
        # Usnap.setUp()
        
        for i in range(self.N):
            lfine.append(self.fineSnapShotList[i].to_petsc().vec())
            lcoarse.append(CoarseSnaps[i].to_petsc().vec())

        CM = self.l2ScalarProductMatrix.to_petsc().mat()

        Udof = lfine[0].copy() # To get the same subdivision in // 

        for i in range(self.N):
            CM.mult(lfine[i],Udof)
            self.reducedBasis.multTranspose(Udof,Usnap)  # ??
            Bh[i,:] = Usnap.getArray()

            CM.mult(lcoarse[i],Udof)
            self.reducedBasis.multTranspose(Udof,Usnap)  # ??
            BH[i,:] = Usnap.getArray()


        #regularization (AT@A +lambda I_d)^-1
        for i in range(self.N):
            R[i,:]=(np.linalg.inv(BH.transpose()@BH+lambd*np.eye(self.N))@BH.transpose()@Bh[:,i])

        return R


    """
    Handle Gram-Schmidt orthogonalization
    """
    def scalarL2(self, u, v):
        """Return the ernegy scalar product associed to the L2 scalar product matrix (mass matrix)
                    int_X(u v)

        Args:
            u (PETSc.Vec): vector
            v (PETSC.Vec): second vector

        Returns:
            float: v.T @ ML2 @ u
        """
        return v.dot( self.l2ScalarProductMatrix.to_petsc().mat() * u )   # v.T @ A @ u

    def scalarH1(self, u, v):
        """Return the ernegy scalar product associed to the H1 scalar product matrix
                    int_X(\nabla u \nabla v)

        Args:
            u (PETSc.Vec): vector
            v (PETSC.Vec): second vector

        Returns:
            float: v.T @ MH1 @ u
        """
        return v.dot( self.h1ScalarProductMatrix.to_petsc().mat() * u )   # v.T @ A @ u

    def normL2(self, u):
        """Compute the L2 norm of the given vector

        Args:
            u (PETSc.Vec): vector

        Returns:
            float: ||u||_L2
        """
        return np.sqrt(self.scalarL2(u, u))

    def normH1(self, u):
        """Compute the H1 norm of the given vector

        Args:
            u (PETSc.Vec): vector

        Returns:
            float: ||u||_H1
        """
        return np.sqrt(self.scalarH1(u, u))

    def orthonormalizeH1(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using H1 norm
        (the optional argument is not needed)
        """
        ub,vb = self.reducedBasis.createVecs()
        Z = []
        for i in range(self.N):
            ub[:] = self.reducedBasis[:,i]
            Z.append(ub)

        # Z[0] /= self.normH1(Z[0])

        for n in range(0, len(Z)):
            s = Z[0].duplicate()
            s.set(0)
            for m in range(n):
                s += self.scalarH1(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp / self.normH1(z_tmp)

        for i in range(self.N):
            self.reducedBasis[:,i] = Z[i][:]

        self.reducedBasis.assemble()

        if not (self.checkH1Orthonormalized() ) and nb < 10:
            self.orthonormalizeH1(nb=nb+1)
        elif rank == 0:
            # pass
            print(f"[NIRB] Gram-Schmidt H1 orthonormalization done after {nb+1} step"+['','s'][nb>0])


    def orthonormalizeL2(self, nb=0):
        """Use Gram-Schmidt algorithm to orthonormalize the reduced basis using L2 norm
        (the optional argument is not needed)
        """
        ub,vb = self.reducedBasis.createVecs()
        Z = []
        for i in range(self.N):
            ub[:] = self.reducedBasis[:,i]
            Z.append(ub)

        # Z[0] /= self.normL2(Z[0])

        for n in range(0, len(Z)):
            s = Z[0].duplicate()
            s.set(0)
            for m in range(n):
                s += self.scalarL2(Z[n], Z[m]) * Z[m]
            z_tmp = Z[n] - s
            Z[n] = z_tmp / self.normL2(z_tmp)

        for i in range(self.N):
            self.reducedBasis[:,i] = Z[i][:]

        self.reducedBasis.assemble()

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
        h1ScalPetsc = self.h1ScalarProductMatrix.to_petsc().mat()
        print('rank ', rank, 'l2 mat size', h1ScalPetsc.local_size, 'mat basis size', self.reducedBasis.local_size)
        matH1 = h1ScalPetsc.PtAP(self.reducedBasis)

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

    def checkL2Orthonormalized(self, tol=1e-8):
        """Check if the reduced basis is L2 orthonormalized.

        Args:
            tol (float, optional): Tolerance. Defaults to 1e-8.

        Returns:
            bool: True if the reduced basis is L2 orthonormalized
        """
        l2ScalPetsc = self.l2ScalarProductMatrix.to_petsc().mat()
        
        print('rank ', rank, 'l2 mat size', l2ScalPetsc.local_size, 'mat basis size', self.reducedBasis.local_size)
        matL2 = l2ScalPetsc.PtAP(self.reducedBasis)

        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    if abs(matL2[i, j] - 1) > tol:
                        return False
                    # assert abs(matL2[i, j] - 1) < tol, f"L2 [{i}, {j}] {matL2[i, j]}"
                    # assert abs(matH1[i, j] - 1) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"
                else:
                    if abs(matL2[i, j]) > tol :
                        return False
                    # assert abs(matL2[i, j]) < tol, f"L2 [{i}, {j}] {matL2[i, j]}"
                    # assert abs(matH1[i, j]) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"
        return True

    def saveData(self, path="./"):
        """Save the data generated by the offline phase

        Args:
            path (str, optional): Path where files are saved. Defaults to "./".
        """
        if not os.path.exists(path):
            os.makedirs(path)
        SavePetscArrayBin(os.path.join(path, "massMatrix.dat"), self.l2ScalarProductMatrix.to_petsc().mat())
        SavePetscArrayBin(os.path.join(path, "stiffnessMatrix.dat"), self.h1ScalarProductMatrix.to_petsc().mat())
        SavePetscArrayBin(os.path.join(path, "reducedBasisU.dat"), self.reducedBasis)
        if self.doRectification:
            SavePetscArrayBin("rectificationMatrix.dat", self.RectificationMat)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Data saved in {os.path.abspath(path)}")



### ONLINE PHASE ###

class nirbOnline(ToolboxModel):

    def __init__(self, dimension, H, h, toolboxType, cfg_path, model_path, geo_path,
        order=1, doRectification=True) -> None:
        """Initialize the NIRB online class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxType (str): toolbox used (heat or fluid)
            cfg_path (str): path to cfg file
            model_path (str): path to json file
            geo_path (str): path to geo file
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """

        super().__init__(dimension, H, h, toolboxType, cfg_path, model_path, geo_path, order)

        self.doRectification = doRectification

        self.l2ScalarProductMatrix = None
        self.h1ScalarProductMatrix = None
        self.reducedBasis = None 
        self.N = 0

        super().initCoarseToolbox()

        self.interpolationOperator = self.createInterpolator(self.tbCoarse, self.tbFine)
        self.exporter = None

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Initialization done")


    # def computeErrors(self, mu=None, exporter=None):
    #     """Compute errors between nirb solution and FE solution computed in fine mesh

    #     Args:
    #         mu (ParameterSpaceElement) : parameter
    #         exporter (feelpp.exporter) : Exporter to export data for visualization

    #     Returns:
    #         list: containing
    #             - the size of the reduced basis N
    #             - the errors Linf and L2 between nirb solution and FE fine solution
    #             - the errors Linf and L2 between inteprolated FE coarse solution and nirb solution
    #             - the errors Linf and L2 between inteprolated FE coarse solution and FE fine solution
    #     """
    #     if mu == None:
    #         mu =self.onlineParam

    #     fineSol = self.getToolboxSolution(self.tbFine, mu)
    #     # fineSol = self.interpSol

    #     error = []
    #     error.append(self.N)

    #     # norm(U_h - U_h^N)
    #     diffSol = (fineSol - self.onlineSol).to_petsc().vec()
    #     error.append(diffSol.norm())
    #     error.append(diffSol.norm(PETSc.NormType.NORM_INFINITY))

    #     # norm(U_H - U_h)
    #     diffSol = (fineSol - self.interpSol).to_petsc().vec()
    #     error.append(diffSol.norm())
    #     error.append(diffSol.norm(PETSc.NormType.NORM_INFINITY))

    #     # norm(U_H - U_h^N)
    #     diffSol = (self.interpSol - self.onlineSol).to_petsc().vec()
    #     error.append(diffSol.norm())
    #     error.append(diffSol.norm(PETSc.NormType.NORM_INFINITY))


        # # # Export fine solution
        # if exporter is not None:
        #     if self.order==1:
        #         exporter.addP1c("U_fine", fineSol)
        #     elif self.order==2:
        #         exporter.addP2c("U_fine", fineSol)

        # return error

    def getCompressedSol(self, mu=None, solution=None):
        """
        get the projection of given solution from fine mesh into the reduced space

        Parameters
        ----------
        solution (feelpp._discr.Element) : the solution to be projected
        mu (ParameterSpaceElement) : parameter to compute the solution if not given 

        return :
        compressedSol (petsc.Vec) : the compressed solution, of size (self.N)
        """
        assert (mu != None) or (solution != None), f"One of the arguments must be given: solution or mu"

        if solution is None:
            sol = self.getInterpSol(mu)
        else :
            sol = solution

        compressedSol,ur = self.reducedBasis.createVecs()      # Get vectors with the same parallel layout as the matrix

        self.l2ScalarProductMatrix.mult(sol.to_petsc().vec(),ur)
        self.reducedBasis.multTranspose(ur, compressedSol)           # ???

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
        resPETSc = self.Xh.element().to_petsc()

        # resPETSc, u = self.reducedBasis.createVecs()  # ??

        compressedSol = self.getCompressedSol(mu)

        if self.doRectification:
            coef = compressedSol.copy()
            self.RectificationMat.mult(compressedSol,coef)
            self.reducedBasis.mult(coef, resPETSc.vec())
            print("[NIRB] Solution computed with Rectification post-process ")
        else :
            self.reducedBasis.mult(compressedSol, resPETSc.vec())

        onlineSol = self.Xh.element(resPETSc)

        # to export, call self.exportField(onlineSol, "U_nirb")

        # No raison to return resPETSc but it makes it work (TO CHECK)
        return onlineSol, resPETSc

    def loadData(self, path="./"):
        """Load the data generated by the offline phase

        Args:
            path (str, optional): Path where files are saved. Defaults to "./".
        """
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Loading data from {os.path.abspath(path)}")
        self.l2ScalarProductMatrix = LoadPetscArrayBin(os.path.join(path, "massMatrix.dat"))
        self.l2ScalarProductMatrix.assemble()
        self.h1ScalarProductMatrix = LoadPetscArrayBin(os.path.join(path, "stiffnessMatrix.dat"))
        self.h1ScalarProductMatrix.assemble()
        self.reducedBasis = LoadPetscArrayBin(os.path.join(path, "reducedBasisU.dat"))
        self.reducedBasis.assemble()
        self.N = self.reducedBasis.size[1]
        if self.doRectification:
            self.RectificationMat = LoadPetscArrayBin("rectificationMatrix.dat")
            self.RectificationMat.assemble()

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
            print("Exporter not initialized, pease call initExporter() first")

    def saveExporter(self):
        """save the exporter to disk
        """
        if self.exporter is not None:
            self.exporter.save()
        else:
            print("Exporter not initialized, please call initExporter() first")