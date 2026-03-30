import feelpp.core as fppc
import pytest
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heat import *
from _heat_test_utils import run_heat_case
from _case_paths import toolbox_case

heat_cases = [
    (toolbox_case('heat/Building/ThermalBridgesENISO10211/case2.cfg'), 2, 1),
    (toolbox_case('heat/Building/ThermalBridgesENISO10211/case2.cfg'), 2, 2),
]


@pytest.mark.parametrize("casefile,dim,order", heat_cases)
def test_heat(casefile,dim,order):
    run_heat_case(casefile, dim, order)

parts = [2,3,6]


@pytest.mark.skip(reason="no way of currently testing this")
@pytest.mark.parametrize("nparts", parts)
def test_heat_ensemble(nparts):
    casefile = toolbox_case('heat/Building/ThermalBridgesENISO10211/case2.cfg')
    c=fppc.readCfg(casefile)
    fppc.Environment.setConfigFile(casefile)
    dim = int(c['feelpp']['case.dimension'])
    order = [int(s) for s in c['feelpp']['case.discretization'].split("P") if s.isdigit()][0]
    if fppc.Environment.numberOfProcessors() > 1 and fppc.Environment.numberOfProcessors() % nparts == 0:
        color, w, wglob = fppc.Environment.worldCommPtr().split(nparts)
        fppc.Environment.changeRepository(c['feelpp']['directory']+"/{}_{}".format(nparts,color))
        f = heat(dim=dim, order=order, worldComm=w)
        if not f.isStationary():
            f.setTimeFinal(10*f.timeStep())
        #simulate(f)
        # meas = f.postProcessMeasures().values()
        # try:
        #     import pandas as pd
        # 
        #     df = pd.DataFrame([meas])
        # except ImportError:
        #     print("cannot import pandas, no problem it was just a test")
        # 
        # return not f.checkResults()

def test_heat_alg():
    fppc.Environment.setConfigFile(
        toolbox_case('heat/Building/ThermalBridgesENISO10211/case2.cfg'))
    f = heat(dim=2, order=1)
    f.init()
    if f.isStationary():
        f.solve()
        maf=f.algebraicFactory()
        A=maf.matrix().to_petsc().mat()
        b=maf.rhs().to_petsc().vec()
        from petsc4py import PETSc
        KSP_TYPE = PETSc.KSP.Type.GMRES
        PC_TYPE = PETSc.PC.Type.LU

        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_SELF)
        ksp.setType(KSP_TYPE)
        reshist = {}

        def monitor(ksp, its, rnorm):
            reshist[its] = rnorm
            #print("[petsc4py] Iteration {} Residual norm: {}".format(its,rnorm))
        ksp.setMonitor(monitor)
        pc = ksp.getPC()
        pc.setType(PC_TYPE)
        ksp.setOperators(A)
        ksp.setConvergenceHistory()
        x=b.duplicate()
        x.set(0)
        ksp.solve(b, x)

        
        T=f.fieldTemperature().to_petsc().vec()
        print("max norm T:", T.norm(PETSc.NormType.NORM_INFINITY))
        err=T-x
        n2 = err.norm(PETSc.NormType.NORM_2)
        print("error norm:",n2)
        assert(n2<1e-8)


