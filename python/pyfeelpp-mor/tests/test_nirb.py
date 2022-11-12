import sys
import os
import pytest


import feelpp
from feelpp.mor.nirb.nirb import *

# desc : ((toolboxtype, 'model_directory', cfg, json, geo, H, h, dimension, doRectification), 'name-of-the-test')
cases = [
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', False), 'lid-driven-cavity w/o rect.'),
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', True) , 'lid-driven-cavity rect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', False, False), 'square2d w/o rect wogreedy'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', True, False) , 'square2d rect wogreedy'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', True, True) , 'square2d rect egreedy'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', False, False), 'thermal-fin-3d w/o rect wogreedy'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', True, False) , 'thermal-fin-3d rect wogreedy'),
        ]
cases_params, cases_ids = list(zip(*cases))


def run_offline(model_path, rect, greedy):
    nbSnap = 6
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect
    nirb_off = nirbOffline(**nirb_config, initCoarse=True)
    nirb_off.generateOperators(coarse=True)

    if greedy:
        _ = nirb_off.initProblemGreedy(100, 1e-5, Nmax=nbSnap, computeCoarse=True, samplingMode="random")
    else:
        nirb_off.initProblem(nbSnap)
    nirb_off.generateReducedBasis(regulParam=1.e-10)

    nirb_off.saveData(force=True)

    assert nirb_off.checkL2Orthonormalized(), "L2 orthonormalization failed"
    # assert nirb_off.checkH1Orthonormalized(), "H1 orthonormalization failed"


def run_online(model_path, rect):
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect
    nirb_on = nirbOnline(**nirb_config)
    err = nirb_on.loadData()
    assert err == 0, "loadData failed"

    mu = nirb_on.Dmu.element()

    uHh = nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)


@pytest.mark.parametrize("dir,cfg,json,rect,greedy", cases_params, ids=cases_ids)
def test_nirb(dir, cfg, json, rect, greedy, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    run_offline(model_path, rect, greedy)
    run_online(model_path, rect)