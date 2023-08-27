import sys
import os
import pytest

import feelpp
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.greedy import *
from feelpp.mor.nirb.nirbOffline import run_offline_pod, run_offline_greedy
from feelpp.mor.nirb.nirbOnline import run_online

# desc : (('path', 'config-file', 'model-file', rectification), 'name-of-the-test')
casesNirb = [
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', False), 'lid-driven-cavity w/o rect.'),
        #  (('testcase/nirb/lid-driven-cavity/', 'cfd2d.cfg', 'cfd2d.json', True) , 'lid-driven-cavity rect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', False) , 'square2d drect'),
         (('testcase/nirb/square', 'square.cfg', 'square.json', True) , 'square2d erect'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json', True) , 'thermal-fin-3d rect'),
        ]
# NB: for the name of the test, erect is a keyword standing for "enable rectification", and drect for "disable rectification"
cases_params_nirb, cases_ids_nirb = list(zip(*casesNirb))

casesInit = [
         (('testcase/nirb/square', 'square.cfg', 'square.json') , 'square-2D'),
         (('testcase/nirb/thermal-fin-3d', 'thermal-fin.cfg', 'thermal-fin.json'), 'thermal-fin-3d'),
        ]
cases_paramsInit, cases_idsInit = list(zip(*casesInit))



def check_saved(config_nirb, loadPath, nirb_off):
    print("\n\nCheck saved data")
    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel()
    err = nirb_on.loadData(path=loadPath)
    assert err == 0

    assert nirb_off.N == nirb_on.N, "N is not the same"
    for i in range(nirb_off.N):
        norm_offline = nirb_off.reducedBasis[i].l2Norm()
        norm_online = nirb_on.reducedBasis[i].l2Norm()
        assert abs(norm_offline - norm_online) < 1e-10, f"norms are not the same for i={i} : {norm_offline} != {norm_online}"



@pytest.mark.parametrize("dir, cfg, json, rect", cases_params_nirb, ids=cases_ids_nirb)
def test_nirb_pod(dir, cfg, json, rect, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)
    feelpp.Environment.setConfigFile(casefile)

    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = rect
    nirb_config['doGreedy'] = False
    # It is just for tests so we take large mesh to fasten the computation
    if nirb_config['dim'] == 2:
        nirb_config['H'] = 0.1
        nirb_config['h'] = 0.01
    else:
        nirb_config['H'] = 0.5
        nirb_config['h'] = 0.2
    nirb_config['nbSnapshots'] = 10

    nirb_off = run_offline_pod(nirb_config)
    path = nirb_off.saveData(force=True)

    check_saved(nirb_config, path, nirb_off)

    s = nirb_off.Dmu.sampling()
    s.sample(10, "random")

    Nbasis = nirb_off.N

    # Check that the online solution is indeed computed
    run_online(nirb_config, path, Xi=s.getVector(), Nb=Nbasis, rectification=rect)
    run_online(nirb_config, path, Xi=s.getVector(), rectification=rect)

    if Nbasis > 4:
        # Check that we can load smaller basis
        run_online(nirb_config, path, nbSnap=Nbasis-2, Xi=s.getVector(), rectification=rect)
        # Check that we can compute solution with a subbasis
        run_online(nirb_config, path, Xi=s.getVector(), Nb=Nbasis - 2, rectification=rect)

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
    nirb_config['doGreedy'] = True
    # It is just for tests so we take large mesh to fasten the computation
    if nirb_config['dim'] == 2:
        nirb_config['H'] = 0.1
        nirb_config['h'] = 0.01
    else:
        nirb_config['H'] = 0.5
        nirb_config['h'] = 0.2

    nirb_offline, _, _ = run_offline_greedy(nirb_config, 5, 200, Nmax=20)
    path = nirb_offline.saveData(force=True)
    check_saved(nirb_config, path, nirb_offline)

    Nbasis = nirb_offline.N
    s = nirb_offline.Dmu.sampling()
    s.sample(10, "random")

    # Check that the online solution is indeed computed
    run_online(nirb_config, path, Xi=s.getVector(), Nb=Nbasis, rectification=rect)
    run_online(nirb_config, path, Xi=s.getVector(), rectification=rect)

    if Nbasis > 4:
        # Check that we can load smaller basis
        run_online(nirb_config, path, nbSnap=Nbasis-2, Xi=s.getVector(), rectification=rect)
        # Check that we can compute solution with a subbasis
        run_online(nirb_config, path, Xi=s.getVector(), Nb=Nbasis - 2, rectification=rect)

@pytest.mark.parametrize("dir, cfg, json", cases_paramsInit, ids=cases_idsInit)
def test_initializer(dir, cfg, json, init_feelpp):
    e = init_feelpp
    casefile = os.path.join(os.path.dirname(__file__), dir, cfg)
    model_path = os.path.join(os.path.dirname(__file__), dir, json)

    feelpp.Environment.setConfigFile(casefile)
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = True
    nirb_config['doGreedy'] = True
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
