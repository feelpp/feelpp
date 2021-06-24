import sys
import feelpp
from feelpp.toolboxes.heat import *

def test_heat():
    feelpp.Environment.setConfigFile(
        'heat/Building/ThermalBridgesENISO10211/thermo2dCase2.cfg')
    f = heat(dim=2, order=1)
    f.init()
    #f.printAndSaveInfo()
    if f.isStationary():
        f.solve()
        f.exportResults()
    else:
        if not f.doRestart():
            f.exportResults(f.timeInitial())
        while not f.timeStepBase().isFinished():
            if f.worldComm().isMasterRank():
                print("============================================================\n")
                print("time simulation: ", f.time(), "s \n")
                print("============================================================\n")
            f.solve()
            f.exportResults()
            f.updateTimeStep()
    return not f.checkResults()


def test_heat_alg():
    feelpp.Environment.setConfigFile(
        'heat/Building/ThermalBridgesENISO10211/thermo2dCase2.cfg')
    f = heat(dim=2, order=1)
    f.init()
    if f.isStationary():
        maf=f.algebraicFactory()
        A=maf.matrix().mat()
        b=maf.rhs().vec()
        from petsc4py import PETSc
        KSP_TYPE = PETSc.KSP.Type.PREONLY
        PC_TYPE = PETSc.PC.Type.LU
        
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_SELF)
        ksp.setType(KSP_TYPE)
        reshist = {}

        def monitor(ksp, its, rnorm):
            reshist[its] = rnorm
            print("Iteration {} Residual norm: {}".format(its,rnorm))
        ksp.setMonitor(monitor)
        pc = ksp.getPC()
        pc.setType(PC_TYPE)
        ksp.setOperators(A)
        ksp.setConvergenceHistory()
        x=b.duplicate()

        ksp.solve(b, x)

        f.solve()
        T=f.fieldTemperature().to_petsc()
        err=T-x
        n2 = err.norm(PETSc.NormType.NORM_2)
        assert(n2<1e-10)





