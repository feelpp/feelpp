import sys
import os
import pytest

import feelpp
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.greedy import *

# desc : ((toolboxtype, 'model_directory', cfg, json, geo, H, h, dimension, doRectification), 'name-of-the-test')
casesNirb = [
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', False), 'lid-driven-cavity w/o rect.'),
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', True) , 'lid-driven-cavity rect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', False), 'square2d w/o rect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', True) , 'square2d rect'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', False), 'thermal-fin-3d w/o rect'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', True) , 'thermal-fin-3d rect'),
        ]
# NB: for the name of the test, wogreedy is a keyword standing for "without greedy", and egreedy for "enable greedy"
cases_params_nirb, cases_ids_nirb = list(zip(*casesNirb))

casesInit = [
         (('testcase/nirb/square', 'square.cfg', 'square.json') , 'square-2D'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json'), 'thermal-fin-3d'),
        ]
cases_paramsInit, cases_idsInit = list(zip(*casesInit))


def run_offline(model_path, rect):
    nbSnap = 6
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect
    nirb_off = nirbOffline(**nirb_config, initCoarse=True)
    nirb_off.initModel()
    nirb_off.generateOperators(coarse=True)

    _ = nirb_off.initProblem(nbSnap)
    RIC = nirb_off.generateReducedBasis(regulParam=1.e-10)

    tolortho =1.e-8
    nirb_off.orthonormalizeL2(tol=tolortho)

    assert nirb_off.checkL2Orthonormalized(tol=tolortho), "L2 orthonormalization failed"
    # assert nirb_off.checkH1Orthonormalized(), "H1 orthonormalization failed"

    nirb_off.saveData(force=True)



def run_offline_greedy(model_path, rect):
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect

    tb = ToolboxModel(**nirb_config)
    tb.initModel(initCoarse=True)
    interpolator = tb.createInterpolator(tb.tbCoarse, tb.tbFine)

    nirb_off = nirbOffline(initCoarse=True, **nirb_config)
    nirb_off.setModel(tb)

    nirb_on = nirbOnline(**nirb_config)
    nirb_on.setModel(tb, interpolationOperator=interpolator)

    nirb_off.generateOperators(coarse=True)

    res = initProblemGreedy(nirb_off, nirb_on, 5, 1000, 1e-5, Nmax=10)



def run_online(model_path, rect):
    nbSnap=6
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect
    nirb_on = nirbOnline(**nirb_config)
    nirb_on.initModel()
    err = nirb_on.loadData(nbSnap=nbSnap)
    assert err == 0, "loadData failed"

    mu = nirb_on.Dmu.element()

    uHh = nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)


@pytest.mark.parametrize("dir,cfg,json,rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    run_offline(model_path, rect)
    run_online(model_path, rect)

@pytest.mark.parametrize("dir,cfg,json", cases_paramsInit, ids=cases_idsInit)
def test_nirb_greedy(dir, cfg, json, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    run_offline_greedy(model_path, True)
    run_online(model_path, True)


@pytest.mark.parametrize("dir,cfg,json", cases_paramsInit, ids=cases_idsInit)
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
