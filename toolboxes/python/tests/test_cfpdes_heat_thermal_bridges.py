from _cfpdes_test_utils import run_cfpdes_case


def test_cfpdes_heat_thermal_bridges():
    run_cfpdes_case("cfpdes/heat/ThermalBridgesENISO10211/thermo2dCase2.cfg", 2)
