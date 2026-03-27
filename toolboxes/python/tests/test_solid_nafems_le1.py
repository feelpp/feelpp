from _solid_test_utils import run_solid_case


def test_solid_nafems_le1():
    run_solid_case("solid/NAFEMS-LE1/le1.cfg", 2, 2)
