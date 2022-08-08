import sys, os
from nirb import *


if __name__ == "__main__":

    PWD = os.getcwd()

    # set the feelpp environment
    config = feelpp.globalRepository("nirb")
    e=feelpp.Environment(sys.argv, opts = toolboxes_options("heat"), config=config)

    # fineness of two grids
    H = 0.1
    h = H**2

    # load the model
    cfg_path = os.path.join(PWD, "nirb/square/square.cfg")
    geo_path = os.path.join(PWD, "nirb/square/square.geo")
    model_path = os.path.join(PWD,"nirb/square/square.json")

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

    print("Size of field temperature :")
    print("    fine :", temperatureCorase.size())
    print("  coarse :", temperatureFine.size())

    I1 = createInterpolator(tbFine, tbCoarse)
    temperatureInterpolate = I1.interpolate(temperatureFine)

    print("Interpolate fine to coarse :")
    print("size :", temperatureInterpolate.size())
    print("Norm :", (temperatureInterpolate - temperatureCorase).to_petsc().vec().norm() )

    I2 = createInterpolator(tbCoarse, tbFine)
    temperatureCoraseInterpolate = I2.interpolate(temperatureCorase)

    print("Interpolate coarse to fine :")
    print('size :', temperatureCoraseInterpolate.size())
    print('Norm :', (temperatureCoraseInterpolate - temperatureFine).to_petsc().vec().norm())
