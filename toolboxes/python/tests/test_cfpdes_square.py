from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes_square():
    run_cfpdes_case("cfpdes/square/square2d.cfg", 2)
