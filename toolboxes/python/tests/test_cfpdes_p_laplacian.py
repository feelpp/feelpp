from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes_p_laplacian():
    run_cfpdes_case("cfpdes/p-laplacian/regularized.cfg", 2)
