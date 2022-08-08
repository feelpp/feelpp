import sys, os
import feelpp
import feelpp.mor as mor
from feelpp.toolboxes.heat import *
import feelpp.interpolation as fi
import json5 as json



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

def loadParameterSpace(model_path):

    crb_model_properties = mor.CRBModelProperties("", feelpp.Environment.worldCommPtr(), "")
    crb_model_properties.setup(model_path)
    Dmu = feelpp.mor._mor.ParameterSpace.New(crb_model_properties.parameters(), feelpp.Environment.worldCommPtr())
    return Dmu


def assembleToolbox(tb, mu):

    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    tb.updateParameterValues()



def createInterpolator(tbCoarse, tbFine):
    """Create an interpolator between two toolboxes
    
    Args:
        tbCorase (Toolbox): coarse toolbox
        tbFine (Toolbox): fine toolbox
    """
    Vh_coarse = tbCoarse.spaceTemperature()
    Vh_fine = tbFine.spaceTemperature()
    interpolator = fi.interpolator(domain = Vh_coarse, image = Vh_fine, range = tbCoarse.rangeMeshElements())
    return interpolator




if __name__ == "__main__":

    PWD = os.getcwd()

    # set the feelpp environment
    config = feelpp.globalRepository("nirb")
    e=feelpp.Environment(sys.argv, opts = toolboxes_options("heat"), config=config)

    # fineness of two grids
    H = 0.1
    h = H**2

    # load the model
    cfg_path = f"{PWD}/nirb/square/square.cfg"
    geo_path = f"{PWD}/nirb/square/square.geo"
    model_path = f"{PWD}/nirb/square/square.json"

    e.setConfigFile(cfg_path)

    model = loadModel(model_path)

    tbCoarse = setToolbox(H, geo_path, model)
    tbFine = setToolbox(h, geo_path, model)

    Dmu = loadParameterSpace(model_path)

    mu = Dmu.element()
    if feelpp.Environment.isMasterRank():
        print(f"Running simulation with mu = {mu}")

    assembleToolbox(tbCoarse, mu)
    assembleToolbox(tbFine, mu)

    tbFine.solve()
    tbCoarse.solve()

    temperatureFine = tbFine.fieldTemperature()
    temperatureCorase = tbCoarse.fieldTemperature()

    print(temperatureCorase.size())
    print(temperatureFine.size())

    I1 = createInterpolator(tbFine, tbCoarse)

    temperatureInterpolate = I1.interpolate(temperatureFine)

    print("size :", temperatureInterpolate.size())

    print("Norm :", (temperatureInterpolate - temperatureCorase).to_petsc().vec().norm() )

    I2 = createInterpolator(tbCoarse, tbFine)

    temperatureCoraseInterpolate = I2.interpolate(temperatureCorase)
    print('size :', temperatureCoraseInterpolate.size())
    print('Norm :', (temperatureCoraseInterpolate - temperatureFine).to_petsc().vec().norm())
