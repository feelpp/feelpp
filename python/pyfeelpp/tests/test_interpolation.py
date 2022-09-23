import sys,os
import feelpp
import pytest
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heat import *
import feelpp.interpolation as fi

cases = [
         (('cases/nirb/square/square.cfg', 'cases/nirb/square/square.geo', 'cases/nirb/square/square.json'), 'square-2d'),
        ]
cases_params, cases_ids = list(zip(*cases))


def setToolbox(h, geo_path, model):

    # load meshes
    mesh_ = feelpp.mesh(dim=2, realdim=2)
    mesh = feelpp.load(mesh_, geo_path, h)

    # set mesh and model properties
    tb = heat(dim=2, order=2)
    tb.setMesh(mesh)
    tb.setModelProperties(model)

    tb.init()

    return tb


def createInterpolator(image_tb,domain_tb):
    """Create an interpolator between two toolboxes
    
    Args:
        source (Toolbox): coarse toolbox
        image (Toolbox): fine toolbox
    """
    Vh_image = image_tb.spaceTemperature()
    Vh_domain = domain_tb.spaceTemperature()
    interpolator = fi.interpolator(domain = Vh_domain, image = Vh_image, range = image_tb.rangeMeshElements())
    return interpolator

@pytest.mark.parametrize("cfg_path, geo_path, model_path", cases_params, ids=cases_ids)
def test_interpolate_constant(cfg_path, geo_path, model_path):

    cfg_path = os.path.join(os.path.dirname(__file__), cfg_path)
    geo_path = os.path.join(os.path.dirname(__file__), geo_path)
    model_path = os.path.join(os.path.dirname(__file__), model_path)

    config = feelpp.globalRepository("nirb")
    e=feelpp.Environment(["pyfeelpp-test-nirb"], opts = toolboxes_options("heat"), config=config)

    # fineness of two grids
    H = 0.1
    h = H**2

    # load the model
    e.setConfigFile(cfg_path)

    model = feelpp.readJson(model_path)

    tbCoarse = setToolbox(H, geo_path, model)
    tbFine = setToolbox(h, geo_path, model)

    I_fineToCoarse = createInterpolator(domain_tb=tbFine, image_tb=tbCoarse)
    I_coarseToFine = createInterpolator(domain_tb=tbCoarse, image_tb=tbFine)

    u_coarse = tbCoarse.spaceTemperature().element()
    u_fine = tbFine.spaceTemperature().element()

    def check_interp(I, u_domain, u_image, u_check_min, u_check_max):
        print("min:{}, max:{}".format(u_domain.min(), u_domain.max()))
        assert abs(u_domain.min()-u_check_min) < 1e-12
        assert abs(u_domain.max()-u_check_max) < 1e-12
        u_image = I.interpolate(u_domain)
        assert abs(u_image.min()-u_check_min) < 1e-12
        assert abs(u_image.max()-u_check_max) < 1e-12

    def check_norm_inter(I, u_domain, u_image_truth):
        u_image = I.interpolate(u_domain)
        diff = (u_image - u_image_truth).to_petsc().vec()
        print("||u_truth - u_inter||2 :", diff.norm())
        assert diff.norm() < 1e-12

    u_coarse.on(feelpp.elements(tbCoarse.mesh()), feelpp.expr("1")) 
    check_interp(I_coarseToFine, u_coarse, u_fine, 1, 1)   
    u_fine.on(feelpp.elements(tbFine.mesh()), feelpp.expr("3")) 
    check_interp(I_fineToCoarse, u_fine, u_coarse, 3, 3)

    u_coarse.on(feelpp.elements(tbCoarse.mesh()), feelpp.expr("x+y:x:y"))
    check_interp(I_coarseToFine, u_coarse, u_fine, 0, 2)
    u_fine.on(feelpp.elements(tbFine.mesh()), feelpp.expr("x+y:x:y"))
    check_interp(I_fineToCoarse, u_fine, u_coarse, 0, 2)

    u_coarse.on(feelpp.elements(tbCoarse.mesh()), feelpp.expr("x^2-3*y^2:x:y"))
    check_interp(I_coarseToFine, u_coarse, u_fine, -3, 1)
    u_fine.on(feelpp.elements(tbFine.mesh()), feelpp.expr("-3*x^2+y^2:x:y"))
    check_interp(I_fineToCoarse, u_fine, u_coarse, -3, 1)

    u_coarse.on(feelpp.elements(tbCoarse.mesh()), feelpp.expr("x^2-3*y^2:x:y"))
    u_fine.on(feelpp.elements(tbFine.mesh()), feelpp.expr("x^2-3*y^2:x:y"))
    check_norm_inter(I_coarseToFine, u_coarse, u_fine)
    check_norm_inter(I_fineToCoarse, u_fine, u_coarse)
