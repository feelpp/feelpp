from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes_heat_thermo2d():
    run_cfpdes_case("cfpdes/heat/thermo2d/testsuite_thermo2d.cfg", 2)
