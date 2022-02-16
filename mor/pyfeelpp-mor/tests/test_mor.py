import sys
import pytest
from petsc4py import PETSc


import feelpp
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *

# desc : (('path/to/cfg/file', dimension, {dict of values in the model}), 'name-of-the-test')
cases = [
         (('thermal-fin/2d/thermal-fin.cfg', 2,
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-2d'),
         (('thermal-fin/3d/thermal-fin.cfg', 3,
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))


def init_toolbox(casefile, dim):

    feelpp.Environment.setConfigFile(casefile)
    feelpp.Environment.changeRepository(casefile)

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    return heatBox



@pytest.mark.parametrize("casefile,dim,mubar_th,mumin_th,mumax_th", cases_params, ids=cases_ids)
def test_init_from_ModelPropeties(casefile, dim, mubar_th, mumin_th, mumax_th,  init_feelpp):
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
    heatBox = init_toolbox(casefile, dim)

    modelParameters = heatBox.modelProperties().parameters()
    Dmu = feelpp.mor._mor.ParameterSpace.New(modelParameters, feelpp.Environment.worldCommPtr())

    mubar = Dmu.mubar()
    mumin = Dmu.mumin()
    mumax = Dmu.mumax()

    mu = Dmu.element()
    print("try to access to the names of elements")
    print(mu.parameterNames())
    print(mumin.parameterNames())


    for p in mubar_th:
        assert( mumin.parameterNamed(p) == mumin_th[p] )
        assert( mumax.parameterNamed(p) == mumax_th[p] )
        assert( mubar.parameterNamed(p) == mubar_th[p] )
