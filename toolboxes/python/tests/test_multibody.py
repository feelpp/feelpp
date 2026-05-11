import numpy as np
import pytest

import feelpp.core as fppc
from _case_paths import toolbox_case

try:
    import feelpp.toolboxes.multibody as multibody_module
    multibody = multibody_module.multibody
    can_import_multibody = multibody_module.has_multibody
except ImportError:
    can_import_multibody = False

try:
    import feelpp.toolboxes.fsi as fsi_module
    fsi = fsi_module.fsi
    can_import_fsi = fsi_module.has_fsi
except ImportError:
    can_import_fsi = False


@pytest.mark.skipif(not can_import_multibody, reason="Required feelpp.toolboxes.multibody module cannot be imported")
def test_multibody_factory_registration(init_feelpp):
    for dim, order_geo in ((2, 1), (2, 2), (3, 1), (3, 2)):
        toolbox = multibody(dim=dim, orderGeo=order_geo, keyword=f"multibody-test-{dim}d-g{order_geo}")
        assert toolbox is not None
        assert toolbox.bodyNames() == []

    with pytest.raises(RuntimeError):
        multibody(dim=4, orderGeo=1)


@pytest.mark.skipif(not can_import_multibody, reason="Required feelpp.toolboxes.multibody module cannot be imported")
def test_multibody_body_inspection(init_feelpp):
    cfg_path = toolbox_case("fsi/magneto/magneto.cfg")
    geo_path = toolbox_case("fsi/magneto/swimmer.geo")
    model_path = toolbox_case("fsi/magneto/magneto.json")

    fppc.Environment.changeRepository(directory="pyfeelpp-tests/multibody/magneto")
    fppc.Environment.setConfigFile(cfg_path)

    mesh = fppc.load(fppc.mesh(dim=2, geo=1, realdim=2), geo_path, 0.1)
    toolbox = multibody(dim=2, orderGeo=1)
    toolbox.setMesh(mesh)
    toolbox.setModelProperties(fppc.readJson(model_path))
    toolbox.init()

    body_names = toolbox.bodyNames()
    assert body_names

    body_name = "fsi-wall"
    assert toolbox.hasBody(body_name)
    assert not toolbox.hasBody("__missing_body__")
    with pytest.raises(IndexError):
        toolbox.body("__missing_body__")

    body = toolbox.body(body_name)
    assert body.mass() > 0
    assert np.asarray(body.massCenter()).size == 2
    assert np.asarray(body.rigidTranslation()).size == 2
    assert np.asarray(body.rigidRotationAngles()).size == 1
    assert np.asarray(body.momentOfInertia_bodyFrame()).size == 1
    assert np.asarray(body.momentOfInertia_inertialFrame()).size == 1
    assert body.fieldDisplacement() is not None
    assert not body.hasElasticDisplacement()
    assert not body.hasElasticVelocity()
    assert body.fieldElasticDisplacement() is None
    assert body.fieldElasticVelocity() is None


@pytest.mark.skipif(not can_import_multibody, reason="Required feelpp.toolboxes.multibody module cannot be imported")
@pytest.mark.skipif(not can_import_fsi, reason="Required feelpp.toolboxes.fsi module cannot be imported")
def test_fsi_model_fluid_exposes_multibody(init_feelpp):
    cfg_path = toolbox_case("fsi/magneto/magneto.cfg")

    fppc.Environment.changeRepository(directory="pyfeelpp-tests/multibody/fsi-magneto")
    fppc.Environment.setConfigFile(cfg_path)

    toolbox = fsi(dim=2, orderU=2, orderP=1, orderGeo=1)
    toolbox.init()

    fluid = toolbox.modelFluid()
    assert fluid is not None

    multibody_model = fluid.multibody()
    assert multibody_model is not None

    body_name = "fsi-wall"
    assert body_name in fluid.bodyNames()
    assert fluid.bodyNames() == multibody_model.bodyNames()
    assert fluid.hasBody(body_name)
    assert not fluid.hasBody("__missing_body__")
    with pytest.raises(IndexError):
        fluid.body("__missing_body__")

    body = fluid.body(body_name)
    assert body.mass() > 0
    assert np.asarray(body.massCenter()).size == 2
    assert body.fieldDisplacement() is not None
