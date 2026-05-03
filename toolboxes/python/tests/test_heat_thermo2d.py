import pytest

from _heat_test_utils import run_heat_case
from _case_paths import toolbox_case

heat_thermo2d_cases = [
    (toolbox_case("heat/thermo2d/thermo2d.cfg"), 2, 1),
    (toolbox_case("heat/thermo2d/thermo2d.cfg"), 2, 2),
]


@pytest.mark.parametrize("casefile,dim,order", heat_thermo2d_cases)
def test_heat_thermo2d(casefile, dim, order):
    run_heat_case(casefile, dim, order)
