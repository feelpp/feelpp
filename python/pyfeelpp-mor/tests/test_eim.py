import sys, os
import pytest
import json5 as json

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import feelpp
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *

from petsc4py import PETSc


# desc : (('path/to/cfg/file', name, dim), 'name-of-the-test')
cases = [
         (('testcase/thermal-fin/2d/thermal-fin.cfg', 'thermal-fin-2d', 2), 'thermal-fin-2d'),
         (('testcase/square/2d/testcase2d.cfg', 'square-2d', 2), 'square-2d'),
         (('testcase/thermal-fin/3d/thermal-fin.cfg', 'thermal-fin-3d', 3), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))



class AffineDecomposition():

    def __init__(self, Aq, Fp) -> None:
        self.Aq = Aq
        self.Fp = Fp

    def computeA(self, betaA):
        A = self.Aq[0].duplicate()

        for q in range(len(self.Aq)):
            A += self.Aq[q] * betaA[q]
        A.assemble()

        return A

    def computeF(self, betaF):
        F = self.Fp[0][0].duplicate()
        F.set(0)

        for p in range(len(self.Fp[0])):
            F += self.Fp[0][p] * betaF[p]

        return F

    def computeLK(self, beta, k):
        F = self.Fp[k][0].duplicate()
        F.set(0)

        for p in range(len(self.Fq[k])):
            F += self.Fp[k][p] * beta[p]

        return F


class EimSolver():

    def __init__(self) -> None:
        self.KSP_TYPE = PETSc.KSP.Type.GMRES
        self.PC_TYPE = PETSc.PC.Type.GAMG
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_SELF)
        self.ksp.setType(self.KSP_TYPE)

    def solve(self, mat, rhs):
        pc = self.ksp.getPC()
        pc.setType(self.PC_TYPE)

        self.reshist = {}

        self.ksp.setOperators(mat)
        self.ksp.setConvergenceHistory()
        sol = rhs.duplicate()
        self.ksp.solve(rhs, sol)

        return sol


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
    for i, F in enumerate(Fq):
        FqP.append([])
        for L in F[0]:
            FqP[i].append(L.to_petsc().vec())
    return FqP



def init_toolboxmor(casefile, name, dim):

    feelpp.Environment.setConfigFile(casefile)
    json_path = feelpp.Environment.expand("$cfgdir")+"/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"

    f = open(json_path, 'r')
    j = json.load(f)
    try:
        j.pop('PostProcess')
    except KeyError as e:
        print(f"There was no section {e} in the model")
    f.close()

    crb_model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    crb_model_properties.setup(json_path)

    model_properties = feelpp.ModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    model_properties.setup(json_path)

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    model = toolboxmor(name=name, dim=dim, time_dependent=False)
    model.setFunctionSpaces( Vh=heatBox.spaceTemperature() )


    def assembleDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleRhs()

    def assembleMDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleMatrix()

    model.setAssembleDEIM(fct=assembleDEIM)
    model.setAssembleMDEIM(fct=assembleMDEIM)
    model.initModel()

    heatBoxDEIM = heat(dim=dim,order=1)
    heatBoxDEIM.setModelProperties(j)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)


    heatBoxMDEIM = heat(dim=dim,order=1)
    heatBoxMDEIM.setModelProperties(j)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineMDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    return model, heatBox, assembleDEIM, assembleMDEIM

@pytest.mark.parametrize("casefile,name,dim", cases_params, ids=cases_ids)
def test_eim(casefile, name, dim, init_feelpp):
    e = init_feelpp

    model, heatBox, assembleDEIM, assembleMDEIM = init_toolboxmor(casefile, name, dim)

    [Aq0, Fq0] = model.getAffineDecomposition()
    Aq = convertToPetscMat(Aq0[0])
    Fq = convertToPetscVec(Fq0)

    AD = AffineDecomposition(Aq, Fq)
    ES = EimSolver()


    print("∥A(mu)∥ , ∥F(mu)∥ , ∥u(mu)∥ ")
    for i in range(10):
        mu = model.parameterSpace().element()
        [betaA, betaF] = model.computeBetaQm(mu)

        F_eim = AD.computeF(betaF[0][0])
        F_tb = assembleDEIM(mu).vec()
        diffF = F_eim - F_tb
        nRelF = diffF.norm()/F_tb.norm()
        # assert(nRelF < 1e-12)


        A_eim = AD.computeA(betaA[0])
        A_tb = assembleMDEIM(mu).mat()
        A_tb.assemble()
        diffA = A_eim - A_tb
        nRelA = diffA.norm()/A_tb.norm()
        # assert(nRelA < 1e-12)

        heatBox.solve()
        # u_tb = heatBox.fieldTemperature().to_petsc().vec()
        u_tb = ES.solve(A_tb, F_tb)
        u_eim = ES.solve(A_eim, F_eim)
        diffU = u_eim - u_tb
        nRelU = diffU.norm()/u_tb.norm()

        print(f"{nRelA:.2e}, {nRelF:.2e}, {nRelU:.2e}")