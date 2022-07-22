import sys,os
import feelpp
import pytest
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heat import *
import feelpp.interpolation as fi
import json5 as json

cases = [
         (('cases/nirb/square/square.cfg', 'cases/nirb/square/square.geo', 'cases/nirb/square/square.json'), 'square-2d'),
        ]
cases_params, cases_ids = list(zip(*cases))


def loadModel(model_path):
    """Load the model from given modle path

    Args:
        model_path (str): path to the model file (JSON)

    Returns:
        json: model loaded
    """
    f = open(model_path, "r")
    model = json.load(f)
    f.close()
    return model


def setToolbox(h, geo_path, model):

    # load meshes
    mesh_ = feelpp.mesh(dim=2, realdim=2)
    mesh = feelpp.load(mesh_, geo_path, h)

    # set mesh and model properties
    tb = heat(dim=2, order=1)
    tb.setMesh(mesh)
    tb.setModelProperties(model)

    tb.init()

    return tb


def createInterpolator(source, image):
    """Create an interpolator between two toolboxes
    
    Args:
        source (Toolbox): coarse toolbox
        image (Toolbox): fine toolbox
    """
    Vh_coarse = source.spaceTemperature()
    Vh_fine = image.spaceTemperature()
    interpolator = fi.interpolator(domain = Vh_coarse, image = Vh_fine, range = source.rangeMeshElements())
    return interpolator

@pytest.mark.parametrize("cfg_path, geo_path, model_path", cases_params, ids=cases_ids)
def test_interpolate_constant(cfg_path, geo_path, model_path):

    PWD=os.getcwd()
    cfg_path = os.path.join(PWD, cfg_path)
    geo_path = os.path.join(PWD, geo_path)
    model_path = os.path.join(PWD, model_path)

    config = feelpp.globalRepository("nirb")
    e=feelpp.Environment(["pyfeelpp-test-nirb"], opts = toolboxes_options("heat"), config=config)

    # fineness of two grids
    H = 0.1
    h = H**2

    # load the model
    e.setConfigFile(cfg_path)

    model = loadModel(model_path)

    tbCoarse = setToolbox(H, geo_path, model)
    tbFine = setToolbox(h, geo_path, model)

    I_fineToCoarse = createInterpolator(tbFine, tbCoarse)
    I_coarseToFine = createInterpolator(tbCoarse, tbFine)

    u = tbCoarse.spaceTemperature().element()
    u.on(feelpp.elements(tbFine.mesh()), feelpp.expr("1"))

    print(u.min(), u.max())


    # failing...
    u_fine = I_fineToCoarse.interpolate(u)
