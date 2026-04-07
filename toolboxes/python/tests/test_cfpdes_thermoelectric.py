from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes_thermoelectric():
    run_cfpdes_case("cfpdes/thermoelectric/ElectroMagnets_HL-31_H1/HL-31_H1.cfg", 3)
