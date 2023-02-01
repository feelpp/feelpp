import sys
import os
import pytest

import feelpp
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.greedy import *
from feelpp.mor.nirb.nirbOffline import run_offline, run_offline_greedy

# desc : (('path', 'config-file', 'model-file', rectification), 'name-of-the-test')
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


@pytest.mark.parametrize("dir, cfg, json, rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect

    nirb_offline = run_offline(nirb_config)
    Nbasis = nirb_offline.N
    # run_online(model_path, rect)

@pytest.mark.parametrize("dir, cfg, json, rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb_greedy(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect

    nirb_offline = run_offline_greedy(nirb_config, 5, 500)
    Nbasis = nirb_offline.N
    # run_online(model_path, True)


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
