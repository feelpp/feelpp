import sys
import feelpp
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heat import *

#def test_heat():
#    feelpp.Environment.setConfigFile(
#        'heat/Building/ThermalBridgesENISO10211/thermo2dCase2.cfg')
#    f = heat(dim=2, order=1)
#    simulate(f)
#    return not f.checkResults()
#

#def test_heat_alg():
#    feelpp.Environment.setConfigFile(
#        'heat/Building/ThermalBridgesENISO10211/thermo2dCase2.cfg')
#    f = heat(dim=2, order=1)
#    f.init()
#    if f.isStationary():
#        f.solve()
#        maf=f.algebraicFactory()
#        A=maf.matrix().to_petsc().mat()
#        b=maf.rhs().to_petsc().vec()
#        from petsc4py import PETSc
#        KSP_TYPE = PETSc.KSP.Type.GMRES
#        PC_TYPE = PETSc.PC.Type.LU
#
#        ksp = PETSc.KSP()
#        ksp.create(PETSc.COMM_SELF)
#        ksp.setType(KSP_TYPE)
#        reshist = {}
#
#        def monitor(ksp, its, rnorm):
#            reshist[its] = rnorm
#            #print("[petsc4py] Iteration {} Residual norm: {}".format(its,rnorm))
#        ksp.setMonitor(monitor)
#        pc = ksp.getPC()
#        pc.setType(PC_TYPE)
#        ksp.setOperators(A)
#        ksp.setConvergenceHistory()
#        x=b.duplicate()
#        x.set(0)
#        ksp.solve(b, x)
#
#        
#        T=f.fieldTemperature().to_petsc().vec()
#        print("max norm T:", T.norm(PETSc.NormType.NORM_INFINITY))
#        err=T-x
#        n2 = err.norm(PETSc.NormType.NORM_2)
#        print("error norm:",n2)
#        assert(n2<1e-8)
#
#



