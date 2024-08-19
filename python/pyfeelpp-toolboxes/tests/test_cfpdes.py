import sys
import pytest
import feelpp.core as fppc
import feelpp.core.quality as fppcq
from feelpp.toolboxes.core import *
from feelpp.toolboxes.cfpdes import *
import pandas as pd

cfpde_cases = [ ('fluid','TurekHron','cfd2.cfg', 2),
                ('.', 'p-laplacian', 'regularized.cfg', 2 ),
                ('.', 'square', 'square2d.cfg', 2), 
                ('thermoelectric', 'ElectroMagnets_HL-31_H1', 'HL-31_H1.cfg', 3),
                ('heat', 'ThermalBridgesENISO10211', 'thermo2dCase2.cfg', 2),
                ('heat', 'thermo2d', 'testsuite_thermo2d.cfg', 2),
                ]

@pytest.mark.parametrize("prefix,case,casefile,dim", cfpde_cases)
def test_cfpdes(prefix,case,casefile,dim):
    fppc.Environment.setConfigFile('cfpdes/{}/{}/{}'.format(prefix,case,casefile))
    f = cfpdes(dim=dim)
    if not f.isStationary():
        # the code below is not working yet see #1763
        f.setTimeFinal( f.timeStep()*10)
    simulate(f)
    return not f.checkResults()
    

# def test_cfpde_nlthermoelectric():
#     prefix, case, casefile, dim= ( 'thermoelectric', 'ElectroMagnets_HL-31_H1', 'HL-31_H1.cfg', 3)
#     fppc.Environment.changeRepository(
#         directory="toolboxes/coefficientformpdes/thermoelectric/ElectroMagnets_HL-31_H1")
#     fppc.Environment.setConfigFile('cfpdes/{}/{}/{}'.format(prefix, case, casefile))
#     linear_case_t = cfpdes(dim=dim)
#     simulate(linear_case_t)
#     assert linear_case_t.checkResults()
#     # now run the nonlinear model using the linear solution as initial guess
#     casefile = 'HL-31_H1_nonlinear.cfg'
#     fppc.Environment.setConfigFile('cfpdes/{}/{}/{}'.format(prefix, case, casefile))
#     nl_case_t = cfpdes(dim=dim)
#     simulate(nl_case_t)
#     return not nl_case_t.checkResults()

def test_cfpdes_remesh():
    fppc.Environment.changeRepository(
        directory="pyfeelpptoolboxes-tests/cfpdes/laplace/l-shape")
    fppc.Environment.setConfigFile('cfpdes/laplace/l-shape/l-shape.cfg')
    f=cfpdes(dim=2)
    simulate(f, export=False)

    e = fppc.exporter(mesh=f.mesh(), name="l-shape", geo="change")
    e.step(0.).setMesh(f.mesh())
    f.exportSolutionToStep( e.step(0.) )
    
    hclose=0.05
    hfar=1

    Xh = fppc.functionSpace(mesh=f.mesh())
    metric = fppc.gradedls(Xh,fppc.boundaryfaces(Xh.mesh()),hclose,hfar)
    e.step(0.).add("metric",metric)
    e.step(0.).add("quality", fppcq.etaQ(f.mesh()))
    e.save()

    new_mesh,cpt = fppc.remesh(
        mesh=f.mesh(), metric="gradedls({},{})".format(hclose, hfar),required_elts=[],required_facets=[],parent=None)
    
    fnew = cfpdes(dim=2)
    fnew.setMesh(new_mesh)
    simulate(fnew,export=False)
    e.step(1.).setMesh(new_mesh)
    fnew.exportSolutionToStep(e.step(1.))

    Xh = fppc.functionSpace(mesh=fnew.mesh())
    metric = fppc.gradedls(Xh,fppc.boundaryfaces(Xh.mesh()),hclose,hfar)
    quality = fppcq.etaQ(fnew.mesh())
    e.step(1.).add("metric",metric)
    e.step(1.).add("quality",quality)
    e.save()
    return not fnew.checkResults()
