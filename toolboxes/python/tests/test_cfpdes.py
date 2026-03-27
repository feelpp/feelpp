from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes():
    run_cfpdes_case("cfpdes/fluid/TurekHron/cfd2.cfg", 2)
