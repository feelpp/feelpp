import sys
import pytest
import feelpp
import feelpp.quality as q
from feelpp.toolboxes.core import *
from feelpp.toolboxes.cfpdes import *
import pandas as pd

cfpde_cases = [ ('fluid','TurekHron','cfd2.cfg', 2),
                ('.', 'p-laplacian', 'regularized', 2 ),
                ('.', 'square', 'square2d', 2), 
                ('thermoelectric', 'ElectroMagnets_HL-31_H1', 'HL-31_H1.cfg', 3),
                #('thermoelectric', 'ElectroMagnets_HL-31_H1','HL-31_H1_nonlinear.cfg', 3),
                ]

@pytest.mark.parametrize("prefix,case,casefile,dim", cfpde_cases)
def test_cfpdes_cfd(prefix,case,casefile,dim):
    feelpp.Environment.changeRepository(
        directory="pyfeelpptoolboxes-tests/cfpdes/{}/{}".format(prefix,case))
    feelpp.Environment.setConfigFile('cfpdes/{}/{}/{}'.format(prefix,case,casefile))
    f = cfpdes(dim=dim)
    simulate(f)
    return not f.checkResults()
    

def test_cfpdes_remesh():
    feelpp.Environment.changeRepository(
        directory="pyfeelpptoolboxes-tests/cfpdes/laplace/l-shape")
    feelpp.Environment.setConfigFile('cfpdes/laplace/l-shape/l-shape.cfg')
    f=cfpdes(dim=2)
    simulate(f, export=False)

    e = feelpp.exporter(mesh=f.mesh(), name="l-shape", geo="change")
    e.step(0.).setMesh(f.mesh())
    f.exportSolutionToStep( e.step(0.) )
    
    hclose=0.05
    hfar=1

    Xh = feelpp.functionSpace(mesh=f.mesh())
    metric = feelpp.gradedls(Xh,feelpp.boundaryfaces(Xh.mesh()),hclose,hfar)
    e.step(0.).add("metric",metric)
    e.step(0.).add("quality", q.etaQ(f.mesh()))
    e.save()

    R = feelpp.remesher(mesh=f.mesh())
    R.setMetric(metric)
    new_mesh = R.execute()

    
    fnew = cfpdes(dim=2)
    fnew.setMesh(new_mesh)
    simulate(fnew,export=False)
    e.step(1.).setMesh(new_mesh)
    fnew.exportSolutionToStep(e.step(1.))

    Xh = feelpp.functionSpace(mesh=fnew.mesh())
    metric = feelpp.gradedls(Xh,feelpp.boundaryfaces(Xh.mesh()),hclose,hfar)
    quality = q.etaQ(fnew.mesh())
    e.step(1.).add("metric",metric)
    e.step(1.).add("quality",quality)
    e.save()
    return not fnew.checkResults()
