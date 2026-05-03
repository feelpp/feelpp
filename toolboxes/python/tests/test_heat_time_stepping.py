import pytest

from _heat_test_utils import run_heat_case
from _case_paths import toolbox_case

heat_time_stepping_cases = [
    (toolbox_case("heat/test_time-stepping/test.cfg"), 2, 1),
    (toolbox_case("heat/test_time-stepping/test.cfg"), 2, 2),
]


@pytest.mark.parametrize("casefile,dim,order", heat_time_stepping_cases)
def test_heat_time_stepping(casefile, dim, order):
    run_heat_case(casefile, dim, order)
