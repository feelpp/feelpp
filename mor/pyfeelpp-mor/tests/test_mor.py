import sys
import pytest
from petsc4py import PETSc


import feelpp
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *

cases = [
         (('thermal-fin', '2d', 'thermal-fin.cfg', 2,
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-2d'),
         (('thermal-fin', '3d', 'thermal-fin.cfg', 3,
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01},
            {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))

mubar_th = {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01}
mumin_th = {"k_1": 0.1, "k_2": 0.1, "k_3": 0.1, "k_4": 0.1, "k_0": 1, "Bi": 0.01}
mumax_th = {"k_1": 10, "k_2": 10, "k_3": 10, "k_4": 10, "k_0": 1, "Bi": 1}


def init_toolbox(prefix, case, casefile, dim):

    feelpp.Environment.setConfigFile(f'{prefix}/{case}/{casefile}')
    feelpp.Environment.changeRepository(directory=f'{prefix}/{case}')

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    return heatBox, dim



@pytest.mark.parametrize("prefix,case,casefile,dim,mubar_th,mumin_th,mumax_th", cases_params, ids=cases_ids)
def test_init_environment(prefix, case, casefile, dim, mubar_th, mumin_th, mumax_th,  init_feelpp):
    e = init_feelpp
    heatBox = init_toolbox(prefix, case, casefile, dim)

    # modelProperties = heatBox.modelProperties()
    model = toolboxmor(dim=dim, time_dependent=False)

    Dmu = model.parameterSpace()
    mubar = Dmu.element()


    decomposition = model.getAffineDecomposition()
    Dmu = model.parameterSpace()
    mubar = Dmu.mubar()
    mumin = Dmu.mumin()
    mumax = Dmu.mumax()
    assert len(decomposition) == 2

    print("\n\n\n\n")
    print(mubar)

    for p in mubar_th:
        assert( mumin.parameterNamed(p) == mumin_th[p] )
        assert( mumax.parameterNamed(p) == mumax_th[p] )
        # assert( mubar.parameterNamed(p) == mubar_th[p] )


    # return decomposition

