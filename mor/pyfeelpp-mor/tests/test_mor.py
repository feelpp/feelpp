import sys
import pytest
from petsc4py import PETSc


import feelpp
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *

cases = [
         (('thermal-fin', '2d', 'thermal-fin.cfg', 2, False, False), 'thermal-fin-2d'),
         (('thermal-fin', '2d', 'thermal-fin.cfg', 2, True, False), 'thermal-fin-2d-cached'),
         (('thermal-fin', '3d', 'thermal-fin.cfg', 3, False, False), 'thermal-fin-3d'),
         (('thermal-fin', '3d', 'thermal-fin.cfg', 3, True, False), 'thermal-fin-3d-cached')
        ]
cases_params, cases_ids = list(zip(*cases))


def init_toolbox(prefix, case, casefile, dim, use_cache):

    # if not use_cache:
    #     print('removing cache...')
    #     try:
    #         shutil.rmtree(FEEL_PATH + '/crbdb')
    #         shutil.rmtree(FEEL_PATH + f'{prefix}/{case}')
    #     except FileNotFoundError:
    #         print('You asked to remove cache, but there was none, so no problem !')

    feelpp.Environment.setConfigFile(f'{prefix}/{case}/{casefile}')
    feelpp.Environment.changeRepository(directory=f'{prefix}/{case}')

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    return heatBox, dim



def init_model(prefix, case, casefile, dim, use_cache, time_dependent):
    heatBox, dim = init_toolbox(prefix, case, casefile, dim, use_cache)
    model = toolboxmor(dim=dim, time_dependent=time_dependent)

    Dmu = model.parameterSpace()
    mubar = Dmu.element()

    modelProperties = heatBox.modelProperties()
    param = modelProperties.parameters()
    for p in param:
        assert(p[1].hasMinMax())    # check that parameters can vary
        mubar.setParameterNamed(p[1].name(), p[1].value())

    model.setFunctionSpaces(Vh=heatBox.spaceTemperature())

    def assembleDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleRhs()

    def assembleMDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleMatrix()

    model.setAssembleDEIM(fct=assembleDEIM)


    model.setAssembleDEIM(fct=assembleDEIM)
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    heatBoxDEIM = heat(dim=dim, order=1)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0, mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)

    heatBoxMDEIM = heat(dim=dim, order=1)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineMDEIM(mu):
        for i in range(0, mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    return heatBox, model, time_dependent, mubar


@pytest.mark.parametrize("prefix,case,casefile,dim,use_cache,time_dependent", cases_params, ids=cases_ids)
def test_init_environment(prefix, case, casefile, dim, use_cache, time_dependent):
    heatBox, model, time_dependent, mubar = init_model(prefix, case, casefile, dim, use_cache, time_dependent)

    decomposition = model.getAffineDecomposition()
    assert len(decomposition) == [2,3][time_dependent]


    # return decomposition

