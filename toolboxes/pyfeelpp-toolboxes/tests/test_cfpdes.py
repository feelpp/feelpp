import sys
import feelpp
import feelpp.quality as q
import feelpp.toolboxes.core as tb
import feelpp.toolboxes.cfpdes as cfpdes
import pandas as pd

def test_cfpdes_cfd():
    feelpp.Environment.changeRepository(
        directory="pyfeelpptoolboxes-tests/cfpdes/fluid/turek-hron")
    feelpp.Environment.setConfigFile('cfpdes/fluid/TurekHron/cfd2.cfg')
    f = cfpdes.cfpdes(dim=2)
    tb.simulate(f)
    
#


#def test_cfpdes_remesh():
#    feelpp.Environment.changeRepository(
#        directory="pyfeelpptoolboxes-tests/cfpdes/laplace/l-shape")
#    feelpp.Environment.setConfigFile('cfpdes/laplace/l-shape/l-shape.cfg')
#    f=cfpdes.cfpdes(dim=2)
#    tb.simulate(f, export=False)
#
#    e = feelpp.exporter(mesh=f.mesh(), name="l-shape", geo="change")
#    e.step(0.).setMesh(f.mesh())
#    f.exportSolutionToStep( e.step(0.) )
#    
#    hclose=0.05
#    hfar=1
#
#    Xh = feelpp.functionSpace(mesh=f.mesh())
#    metric = feelpp.gradedls(Xh,feelpp.boundaryfaces(Xh.mesh()),hclose,hfar)
#    e.step(0.).add("metric",metric)
#    e.step(0.).add("quality", q.etaQ(f.mesh()))
#    e.save()
#
#    R = feelpp.remesher(mesh=f.mesh())
#    R.setMetric(metric)
#    new_mesh = R.execute()
#    
#    
#    fnew = cfpdes.cfpdes(dim=2)
#    fnew.setMesh(new_mesh)
#    tb.simulate(fnew,export=False)
#    e.step(1.).setMesh(new_mesh)
#    fnew.exportSolutionToStep(e.step(1.))
#
#    Xh = feelpp.functionSpace(mesh=fnew.mesh())
#    metric = feelpp.gradedls(Xh,feelpp.boundaryfaces(Xh.mesh()),hclose,hfar)
#    quality = q.etaQ(fnew.mesh())
#    e.step(1.).add("metric",metric)
#    e.step(1.).add("quality",quality)
#    e.save()
#    return not fnew.checkResults()
#