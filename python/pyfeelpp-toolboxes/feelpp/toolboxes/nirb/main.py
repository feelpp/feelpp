import sys, os
from nirb import *
import feelpp.toolboxes.core as core


def test_interpolation(tbCoarse, tbFine, function, type_tb):
    if type_tb == "heat":
        u_coarse = tbCoarse.spaceTemperature().element()
        u_fine = tbFine.spaceTemperature().element()
    elif type_tb == "fluid":
        u_coarse = tbCoarse.spaceVelocity().element()
        u_fine = tbFine.spaceVelocity().element()
    else:
        raise ValueError("type_tb must be 'heat' or 'fluid'")

    u_coarse.on(feelpp.elements(tbCoarse.mesh()), feelpp.expr(function))
    u_fine.on(feelpp.elements(tbFine.mesh()), feelpp.expr(function))

    I = createInterpolator(tbCoarse, tbFine, type_tb=type_tb)
    u_fine_inter = I.interpolate(u_coarse)

    diff = (u_fine - u_fine_inter).to_petsc().vec()
    # print(f"f={function}", diff.norm())
    return diff.norm()


def solve_for_mu(mu, tbCoarse, tbFine, type_tb):
    if feelpp.Environment.isMasterRank():
        print(f"Running simulation with mu = {mu}")

    assembleToolbox(tbCoarse, mu)
    assembleToolbox(tbFine, mu)

    tbFine.solve()
    tbCoarse.solve()

    fieldFine = getField(tbFine, type_tb)
    fieldCoarse = getField(tbCoarse, type_tb)

    print("Size of coarse space :", fieldCoarse.size())
    print("Size of fine space   :", fieldFine.size())

    I1 = createInterpolator(tbFine, tbCoarse, type_tb=type_tb)
    temperatureInterpolate = I1.interpolate(fieldFine)

    print("Interpolation Fine -> Coarse")
    print("  Size of interpolated :", temperatureInterpolate.size())
    assert temperatureInterpolate.size() == fieldCoarse.size()
    print("  Norm :", (temperatureInterpolate - fieldCoarse).to_petsc().vec().norm() )

    I2 = createInterpolator(tbCoarse, tbFine, type_tb=type_tb)
    fieldCoarseInterpolate = I2.interpolate(fieldCoarse)

    print("Interpolation Coarse -> Fine")
    print("  Size of interpolated :", fieldCoarseInterpolate.size())
    assert fieldCoarseInterpolate.size() == fieldFine.size()
    print('  Norm :', (fieldCoarseInterpolate - fieldFine).to_petsc().vec().norm())



if __name__ == "__main__":

    argc = len(sys.argv)
    if argc > 1:
        type_tb = sys.argv[1]
    else:
        type_tb = "heat"

    PWD = os.getcwd()

    # set the feelpp environment
    config = feelpp.globalRepository("nirb")
    e = feelpp.Environment(sys.argv, opts = core.toolboxes_options(type_tb), config=config)

    # fineness of two grids
    H = 0.1
    h = H**2

    # load the model
    if type_tb == "heat":
        cfg_path = os.path.join(PWD, "nirb/square/square.cfg")
        geo_path = os.path.join(PWD, "nirb/square/square.geo")
        model_path = os.path.join(PWD,"nirb/square/square.json")
    else:
        cfg_path = os.path.join(PWD, "nirb/lid-driven-cavity/cfd2d.cfg")
        geo_path = os.path.join(PWD, "nirb/lid-driven-cavity/cfd2d.geo")
        model_path = os.path.join(PWD,"nirb/lid-driven-cavity/cfd2d.json")

    e.setConfigFile(cfg_path)

    model = loadModel(model_path)

    tbCoarse = setToolbox(H, geo_path, model, type_tb=type_tb)
    tbFine = setToolbox(h, geo_path, model, type_tb=type_tb)

    Dmu = loadParameterSpace(model_path)

    # mu = Dmu.element()
    # solve_for_mu(mu, tbCoarse, tbFine, type_tb)

    if type_tb == "heat":
        r = test_interpolation(tbCoarse, tbFine, "x+y:x:y", type_tb)
        assert r < 1e-12
        r = test_interpolation(tbCoarse, tbFine, "x^2+y:x:y", type_tb)
        assert r < 1e-12
        r = test_interpolation(tbCoarse, tbFine, "x^3+y:x:y", type_tb)
        print(f"{H},{h},{r}")
    if type_tb == "fluid":
        r = test_interpolation(tbCoarse, tbFine, "{x, y}:x:y", type_tb)
        assert r < 1e-12
        r = test_interpolation(tbCoarse, tbFine, "{x^2+y,1}:x:y", type_tb)
        assert r < 1e-12
        r = test_interpolation(tbCoarse, tbFine, "{x^3+y, y^3}:x:y", type_tb)
        print(f"{H},{h},{r}")
