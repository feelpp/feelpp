import sys
from pathlib import Path
import pytest
import feelpp.core as fppc
from feelpp.toolboxes.core import *
from feelpp.toolboxes.thermoelectric import *
from _case_paths import toolbox_case

thermoelectric_cases = [
    (toolbox_case('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg'), 3, 1, 1),
    (toolbox_case('thermoelectric/ElectroMagnets/HL-31_H1/HL-31_H1.cfg'), 3, 2, 1),
    *[(toolbox_case('thermoelectric/test/test-thermoelectric.cfg'), 2, order, geo)
      for order, geo in [(1, 3), (2, 1), (2, 3), (3, 2), (4, 1)]],
]


@pytest.mark.parametrize("casefile,dim,order,geo_order", thermoelectric_cases)
def test_thermoelectric(casefile, dim, order, geo_order, case_repository):
    fppc.Environment.setConfigFile(casefile)
    fppc.Environment.changeRepository(directory=case_repository)
    f = thermoelectric(dim=dim, orderPotential=order, orderGeometry=geo_order)
    model_file = str(Path(casefile).with_suffix(".json"))
    if dim == 3 and order > 1:
        # The historical end-flux references reproduce the coarse P1 result.
        # Higher orders check the prescribed natural zero flux on V0/V1;
        # retain the original tolerances and every other reference value.
        properties = fppc.readJson(model_file)
        checkers = properties["PostProcess"]["heat"]["Checkers"]["Measures"]
        for end in ("V0", "V1"):
            checkers[f"Statistics_Flux{end}_integrate"]["value"] = 0
        f.setModelProperties(properties)
    else:
        f.setModelProperties(model_file)
    status, _ = simulate(f)
    assert status
    assert f.mesh() is not None
    assert f.mesh().order() == geo_order
    assert f.heatModel().spaceTemperature().order() == order
    assert f.electricModel().spaceElectricPotential().order() == order
