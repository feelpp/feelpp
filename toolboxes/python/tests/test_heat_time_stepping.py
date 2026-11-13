import pytest

from _heat_test_utils import run_heat_case
from _case_paths import toolbox_case

heat_time_stepping_cases = [
    (toolbox_case("heat/test_time-stepping/test.cfg"), 2, order, geo)
    for order, geo in [(1, 1), (1, 3), (2, 1), (2, 3), (3, 2), (4, 1)]
]


@pytest.mark.parametrize("casefile,dim,order,geo_order", heat_time_stepping_cases)
def test_heat_time_stepping(casefile, dim, order, geo_order, case_repository):
    run_heat_case(casefile, dim, order, geo_order, repository=case_repository)
