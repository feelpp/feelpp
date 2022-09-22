import sys
import os
import pytest


import feelpp
from feelpp.mor.nirb.nirb import *

# desc : ((toolboxtype, 'model_directory', cfg, json, geo, H, h, dimension, doRectification), 'name-of-the-test')
cases = [
        #  (('fluid', 'testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', 'cfd2d.geo', 0.1, 0.5**2, 2, False), 'lid-driven-cavity w/o rect.'),
        #  (('fluid', 'testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', 'cfd2d.geo', 0.1, 0.5**2, 2, True), 'lid-driven-cavity rect'),
         (('heat', 'testcase/nirb/square', 'square.cfg', 'square.json', 'square.geo', 0.1, 0.1**2, 2, False), 'square2d w/o rect'),
         (('heat', 'testcase/nirb/square', 'square.cfg', 'square.json', 'square.geo', 0.1, 0.1**2, 2, True), 'square2d rect'),
        #  (('heat', 'testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', 'fin.geo', 1, 0.05, 3, False), 'thermal-fin-3d w/o rect'),
        #  (('heat', 'testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', 'fin.geo', 1, 0.05, 3, True), 'thermal-fin-3d rect'),
        ]
cases_params, cases_ids = list(zip(*cases))


def run_offline(toolboxType, casefile, model_path, geo_path, dim, H, h, rect):
    nbSnap = 6
    nirb_off = nirbOffline(dim, H, h, toolboxType, casefile, model_path, geo_path, doRectification=rect)

    nirb_off.initProblem(nbSnap)
    nirb_off.generateOperators()
    nirb_off.generateReducedBasis(regulParam=1.e-10)

    nirb_off.saveData()

    assert nirb_off.checkL2Orthonormalized(), "L2 orthonormalization failed"
    # assert nirb_off.checkH1Orthonormalized(), "H1 orthonormalization failed"


def run_online(toolboxType, casefile, model_path, geo_path, dim, H, h, rect):
    nirb_on = nirbOnline(dim, H, h, toolboxType, casefile, model_path, geo_path, doRectification=rect)
    nirb_on.loadData()

    mu = nirb_on.Dmu.element()

    uHh, _ = nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)


@pytest.mark.parametrize("toolboxType,dir,cfg,json,geo,H,h,dim, rect", cases_params, ids=cases_ids)
def test_nirb(toolboxType, dir, cfg, json, geo, dim, H, h, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    geo_path = os.path.join(os.path.dirname(__file__), dir, geo)
    feelpp.Environment.setConfigFile(casefile)

    run_offline(toolboxType, casefile, model_path, geo_path, dim, H, h, rect)
    run_online(toolboxType, casefile, model_path, geo_path, dim, H, h, rect)