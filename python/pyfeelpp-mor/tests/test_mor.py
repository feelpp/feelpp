import sys
import os
import pytest
from petsc4py import PETSc


import feelpp
from feelpp.mor import *

# desc : (('path/to/cfg/file', dimension, {mumin}, {mumax}), 'name-of-the-test')
cases = [
         (('testcase/thermal-fin/2d/thermal-fin.cfg', 2,
           {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
           {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-2d'),
         (('testcase/thermal-fin/3d/thermal-fin.cfg', 3,
           {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
           {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-3d'),
         (('testcase/testcase/test.cfg', 2,
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1},
            {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10}), 'test-case'),
        ]
cases_params, cases_ids = list(zip(*cases))


@pytest.mark.parametrize("casefile,dim,mumin_th,mumax_th", cases_params, ids=cases_ids)
def test_init_from_ModelPropeties(casefile, dim, mumin_th, mumax_th,  init_feelpp):
    """Tests the initialisation from ModelProperties

    Args:
        casefile (string): path to config file
        dim (int): dimention of the case
        mubar_th (dict): dict containing the true values of mu_bar (="value" field in the json)
        mumin_th (dict): true minimal values for parameters (="min" field in the json)
        mumax_th (dict): true maximal values for parameters (="max" fiels in the json)
        init_feelpp (_fixture_):
    """
    e = init_feelpp
    feelpp.Environment.setConfigFile(casefile)

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"
    model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    model_properties.setup(model_path)
    modelParameters = model_properties.parameters()
    Dmu = feelpp.mor._mor.ParameterSpace.New(modelParameters, feelpp.Environment.worldCommPtr())

    # mubar = Dmu.mubar()       # /!\ doesn't exist anymore !
    mumin = Dmu.mumin()
    mumax = Dmu.mumax()

    mu = Dmu.element()
    print("try to access to the names of elements")
    print(mu.parameterNames())
    print(mumin.parameterNames())


    for p in mumin_th:
        assert( mumin.parameterNamed(p) == mumin_th[p] )
        assert( mumax.parameterNamed(p) == mumax_th[p] )


@pytest.mark.parametrize("casefile", [cases_params[-1][0]], ids=[cases_ids[-1]])
def test_sampling(casefile, init_feelpp):
    e = init_feelpp
    feelpp.Environment.setConfigFile(casefile)

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"
    model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    model_properties.setup(model_path)
    modelParameters = model_properties.parameters()
    Dmu = feelpp.mor._mor.ParameterSpace.New(modelParameters, feelpp.Environment.worldCommPtr())

    print(feelpp.mor.__file__)
    Nsamples = 10
    s = Dmu.sampling()
    s.sampling(Nsamples, samplingMode="random")

    assert len(s) == Nsamples , "wrong number of samples"
    v = s.getVector()
    assert len(v) == Nsamples , "wrong number of samples"

    s.writeOnFile("samplingTest.sample")
    s_read  = Dmu.sampling()
    N_read = s_read.readFromFile("samplingTest.sample")
    v_read  = s_read.getVector()

    assert N_read == Nsamples, "wrong number of samples read"

    for i in range(Nsamples):
        for p in Dmu.parameterNames():
            diff = abs(v[i].parameterNamed(p) - v_read [i].parameterNamed(p))
            assert diff < 1e-10 , f"wrong sample read at index {i} : {v[i].parameterNamed(p)} != {v_read[i].parameterNamed(p)}"


@pytest.mark.parametrize("casefile,dim,Mu", [cases_params[-1][:3]], ids=[cases_ids[-1]])
def test_param_not_in_range(casefile, dim, Mu,  init_feelpp):
    e = init_feelpp
    feelpp.Environment.setConfigFile(casefile)

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"
    model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    model_properties.setup(model_path)
    modelParameters = model_properties.parameters()
    Dmu = feelpp.mor._mor.ParameterSpace.New(modelParameters, feelpp.Environment.worldCommPtr())

    mu = Dmu.element()
    # just set values
    mu.setParameters(Mu)

    # check that exception is raised when setting a parameter out of range
    try:
        mu.setParameter(0, 5000)
    except ValueError as e:
        assert "out of range" in str(e)

    try:
        mu.setParameterNamed('Bi', -666)
    except ValueError as e:
        assert "out of range" in str(e)
