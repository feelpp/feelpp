from _fsi_test_utils import run_fsi_case


def test_fsi_turek_hron_fsi1():
    run_fsi_case("fsi/TurekHron/fsi1.cfg", 2, 2, 1, 1)
