import sys
import os
import pytest

import feelpp
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.greedy import *
from feelpp.mor.nirb.nirbOffline import run_offline, run_offline_greedy
from feelpp.mor.nirb.nirbOnline import run_online

# desc : (('path', 'config-file', 'model-file', rectification), 'name-of-the-test')
casesNirb = [
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', False), 'lid-driven-cavity w/o rect.'),
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', True) , 'lid-driven-cavity rect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', False, False), 'square2d w/o rect wogreedy'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', True, False) , 'square2d rect wogreedy'),
        #  (('testcase/nirb/square', 'square.cfg', 'square.json', True, True) , 'square2d rect egreedy'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', False, False), 'thermal-fin-3d w/o rect wogreedy'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', True, False) , 'thermal-fin-3d rect wogreedy'),
        ]
# NB: for the name of the test, wogreedy is a keyword standing for "without greedy", and egreedy for "enable greedy"
cases_params_nirb, cases_ids_nirb = list(zip(*casesNirb))

casesInit = [
         (('testcase/nirb/square', 'square.cfg', 'square.json') , 'square-2D'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json'), 'thermal-fin-3d'),
        ]
cases_paramsInit, cases_idsInit = list(zip(*casesInit))




@pytest.mark.parametrize("dir, cfg, json, rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect

    uHh = nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    errorNirb = nirb_on.normMat(uHh - uh)
    errorInterp = nirb_on.normMat(uH - uh)

    # assert errorNirb<0.08, f"higher nirb error value"
    # assert errorInterp<0.05, f"higher interp error value"

@pytest.mark.parametrize("dir, cfg, json, rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb_greedy(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect

    nirb_offline, _, _ = run_offline_greedy(nirb_config, 5, 200, Nmax=20)
    nirb_offline.saveData(force=True)
    Nbasis = nirb_offline.N
    s = nirb_offline.Dmu.sampling()
    s.sampling(10, "random")
    
    # Check that the online solution is indeed computed
    run_online(nirb_config, nirb_offline.outdir, Nbasis, s.getVector())

    if Nbasis > 4:
        # Check that we can load smaller basis
        run_online(nirb_config, nirb_offline.outdir, Nbasis - 2, s.getVector())
        # Check that we can compute solution with a subbasis
        run_online(nirb_config, nirb_offline.outdir, Nbasis, s.getVector(), Nb=Nbasis - 2)

@pytest.mark.parametrize("dir, cfg, json", cases_paramsInit, ids=cases_idsInit)
def test_initializer(dir, cfg, json, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)

    feelpp.Environment.setConfigFile(casefile)
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = True
    tbModel = ToolboxModel(**nirb_config)
    tbModel.initModel()

    nirb_off = nirbOffline(**nirb_config, initCoarse=True)
    nirb_off.setModel(tbModel)
    assert nirb_off.Xh == tbModel.Xh, "Xh not equal"
    assert nirb_off.tbFine == tbModel.tbFine, "tbFine not equal"

    nirb_on = nirbOnline(**nirb_config)
    nirb_on.setModel(tbModel)
    assert nirb_on.Xh == tbModel.Xh, "Xh not equal"
    assert nirb_on.tbFine == tbModel.tbFine, "tbFine not equal"
