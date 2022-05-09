import os, sys
import pytest
import feelpp.mor.generate_basis as g


#        (( prefix, case, casefile, dim, time_dependant), name     )
cases = [
         (('testcase/thermal-fin', '2d', 'thermal-fin.cfg', 2, False), 'thermal-fin-2d'),
         (('testcase/thermal-fin', '3d', 'thermal-fin.cfg', 3, False), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))

CWD = os.getcwd()



@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_sample(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp

    g.dim = dim
    g.time_dependant = time_dependent
    g.algo = 0
    g.size = 40
    g.case = CWD + "/" + prefix + "/" + case
    g.casefile = casefile
    g.config_file = g.case + "/" + casefile
    g.odir = "sample"

    # compute and save the basis
    g.generate_basis()



@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_greedy(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp
    g.dim = dim
    g.time_dependant = time_dependent
    g.algo = 1
    g.size = 60
    g.case = CWD + "/" + prefix + "/" + case
    g.casefile = casefile
    g.config_file = g.case + "/" + casefile
    g.odir = "greedy"
    g.tol = 1e-3

    # compute and save the basis
    g.generate_basis()


@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_POD(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp

    g.dim = dim
    g.time_dependant = time_dependent
    g.algo = 2
    g.size = 15
    g.case = CWD + "/" + prefix + "/" + case
    g.casefile = casefile
    g.config_file = g.case + "/" + casefile
    g.odir = "/POD"
    g.tol = 1e-3

    # compute and save the basis
    g.generate_basis()