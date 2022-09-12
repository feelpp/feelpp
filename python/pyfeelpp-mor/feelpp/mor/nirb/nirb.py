# -*- coding: utf-8 -*-
## Non intrusive reduced basis method with Feel++
## Thomas Saigre, Ali Elarif
## 09/2022


import feelpp
import feelpp.toolboxes.core as core
from utils import *
import feelpp.mor as mor
import feelpp.operators as FppOp
# from NIRBinitCase import *
import Mordicus.Modules.Cemosis.DataCompressors.SnapshotReducedBasis as SRB


import sys, os

class NIRB():

    def __init__(self, dimension, H, h, toolboxOption, cfg_path, model_path, geo_path,
        method="POD", order=1, doRectification=True) -> None:
        """Initialize the NIRB class

        Args:
            dimension (int): dimension of the case
            H (float): coarse mesh size
            h (float): fine mesh size
            toolboxOption (str): toolbox used (heat or fluid)
            cfg_path (str): path to cfg file
            model_path (str): path to json file
            geo_path (str): path to geo file
            method (str, optional): method used to generate the basis. Defaults to "POD".
            order (int, optional): order of discretization. Defaults to 1.
            doRectification (bool, optional): set rectification. Defaults to True.
        """

        self.e = None      # Feel++ environment

        assert dimension in [2,3]
        self.dimension = dimension
        assert method in ["POD", "Greedy"]
        self.H = H
        self.h = h if h is not None else H**2
        self.method = method
        self.order = order

        self.doRectification = doRectification

        self.toolboxOption = toolboxOption
        self.cfg_path = cfg_path
        self.model_path = model_path
        self.geo_path = geo_path

        self.l2ScalarProductMatrix = None
        self.h1ScalarProductMatrix = None

        self.initFeelpp()
        if self.doRectification:
            self.initRectification()
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Initialization done")


    def initFeelpp(self):
        """Initialize Feel++ environment
        """
        config = feelpp.globalRepository(f"nirb/{self.toolboxOption}")
        opts = core.toolboxes_options(self.toolboxOption).add(mor.makeToolboxMorOptions())
        self.e = feelpp.Environment(sys.argv, config=config, opts=opts)
        self.e.setConfigFile(self.cfg_path)

        self.model = loadModel(self.model_path)
        self.tbFine = setToolbox(self.h, self.geo_path, self.model, dim=self.dimension,
                                      order=self.order, type_tb=self.toolboxOption)
        self.Dmu = loadParameterSpace(self.model_path)

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of nodes on the fine mesh : {self.tbFine.mesh().numGlobalPoints()}")


    def initRectification(self):
        """Initialize the rectification problem
        """
        self.tbCoarse = setToolbox(self.H, self.geo_path, self.model, dim=self.dimension,
                                        order=self.order, type_tb=self.toolboxOption)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of nodes on the coarse mesh : {self.tbCoarse.mesh().numGlobalPoints()}")


    def initProblem(self, numberOfInitSnapshots, samplingMode="log-random", computeCoarse=False):
        """Initialize the problem

        Args:
            numberOfInitSnapshots (int): number of snapshots to use for the initialization
            samplingMode (str, optional): sampling mode in the parameter space. Defaults to "log-random".
            computeCoarse (bool, optional): compute snapshots for coarse toolbox, used for rectification. Defaults to False.
        """
        self.fineSnapShotList = []
        self.coarseSnapShotList = []

        s = self.Dmu.sampling()
        s.sampling(numberOfInitSnapshots, samplingMode)
        vector_mu = s.getVector()

        if computeCoarse:
            for mu in vector_mu:
                if feelpp.Environment.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(getSolution(self.tbFine, mu))
                self.coarseSnapShotList.append(getSolution(self.tbCoarse, mu))

        else:
            for mu in vector_mu:
                if feelpp.Environment.isMasterRank():
                    print(f"Running simulation with mu = {mu}")
                self.fineSnapShotList.append(getSolution(self.tbFine, mu, type_tb=self.toolboxOption))

        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of snapshot computed : {len(self.fineSnapShotList)}" )


    def generateOperators(self):
        """Assemble L2 and H1 operators, stored in self.l2ScalarProduct and self.h1ScalarProduct
        """
        if self.l2ScalarProductMatrix is None or self.h1ScalarProductMatrix is None:
            Vh = feelpp.functionSpace(mesh=self.tbFine.mesh())
            self.l2ScalarProductMatrix = FppOp.mass(test=Vh, trial=Vh, range=feelpp.elements(self.tbFine.mesh()))
            self.h1ScalarProductMatrix = FppOp.stiffness(test=Vh, trial=Vh, range=feelpp.elements(self.tbFine.mesh()))

    def generateReducedBasis(self):
        """Generate the reduced basis, and store it in self.reducedOrderBasisU
        """
        # @ali : in the first code, there was a if calling the same function
        self.reducedOrderBasisU = SRB.PODReducedBasisPETSc(self.fineSnapShotList, self.l2ScalarProductMatrix)
        self.N = self.reducedOrderBasisU.size[0]
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Number of modes : {self.N}")
        if self.doRectification:
            interpolateOperator = createInterpolator(self.tbFine, self.tbCoarse, type_tb=self.toolboxOption)
            coarseSnapInterList = []
            for snap in self.coarseSnapShotList:
                coarseSnapInterList.append(interpolateOperator.interpolate(snap))
            self.RectificationMat = SRB.Rectification(coarseSnapInterList, self.fineSnapShotList,
                                        self.l2ScalarProductMatrix, self.h1ScalarProductMatrix)

    def checkOrthonormalized(self, tol=1e-8):
        """Check if the reduced basis is orthonormalized. Only L2 norm is checked.

        Args:
            tol (float, optional): Tolerance. Defaults to 1e-8.
        """
        l2ScalPetsc = self.l2ScalarProductMatrix.to_petsc().mat()
        h1ScalPetsc = self.h1ScalarProductMatrix.to_petsc().mat()
        # with the current version of petsc, we need to transpose the matrix
        self.reducedOrderBasisU.transpose()
        matL2 = l2ScalPetsc.PtAP(self.reducedOrderBasisU)
        matH1 = h1ScalPetsc.PtAP(self.reducedOrderBasisU)
        self.reducedOrderBasisU.transpose()

        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    assert abs(matL2[i, j] - 1) < tol, f"L2 [{i}, {j}] {matL2[i, j]}"
                    # assert abs(matH1[i, j] - 1) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"
                else:
                    assert abs(matL2[i, j]) < tol, f"L2 [{i}, {j}] {matL2[i, j]}"
                    # assert abs(matH1[i, j]) < tol, f"H1 [{i}, {j}] {matH1[i, j]}"


    def saveData(self):
        """Save the data generated by the offline phase
        """
        SavePetscArrayBin("massMatrix.dat", self.l2ScalarProductMatrix.mat())
        SavePetscArrayBin("stiffnessMatrix.dat", self.h1ScalarProductMatrix.mat())
        SavePetscArrayBin("reducedBasisU.dat", self.reducedOrderBasisU)
        if self.doRectification:
            SavePetscArrayBin("rectificationMatrix.dat", self.RectificationMat)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Data saved in {os.getcwd()}")


    ### ONLINE PHASE ###

    def loadData(self):
        """Load the data generated by the offline phase
        """
        self.l2ScalarProductMatrix = LoadPetscArrayBin("massMatrix.dat")
        self.l2ScalarProductMatrix.assemble()
        self.h1ScalarProductMatrix = LoadPetscArrayBin("stiffnessMatrix.dat")
        self.h1ScalarProductMatrix.assemble()
        self.reducedOrderBasisU = LoadPetscArrayBin("reducedBasisU.dat")
        self.reducedOrderBasisU.assemble()
        self.N = self.reducedOrderBasisU.size[0]
        if self.doRectification:
            self.RectificationMat = LoadPetscArrayBin("rectificationMatrix.dat")
            self.RectificationMat.assemble()
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Data loaded from {os.getcwd()}")

    def solveOnline(self, mu, exporter=None):
        """Solve the online problem with the given parameter mu
        """
        coarseSol = getSolution(self.tbCoarse, mu)
        interpolationOperator = createInterpolator(self.tbCoarse, self.tbFine, type_tb=self.toolboxOption)
        interpolatedSol = interpolationOperator.interpolate(coarseSol)

        # Export
        if exporter is not None:
            if self.order == 1:
                exporter.addP1c("U_interp", interpolatedSol)
            elif self.order == 2:
                exporter.addP2c("U_interp", interpolatedSol)
        return interpolatedSol