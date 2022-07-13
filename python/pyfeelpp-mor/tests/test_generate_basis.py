import os, sys
import pytest
import feelpp.mor.generate_basis as g


#        (( prefix, case, casefile, dim, time_dependant), name     )
cases = [
         (('testcase', 'thermal-fin/2d', 'thermal-fin.cfg', 2, False), 'thermal-fin-2d'),
         (('testcase', 'thermal-fin/3d', 'thermal-fin.cfg', 3, False), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))

CWD = os.getcwd()



@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_sample(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp

    dim = dim
    time_dependant = time_dependent
    algo = 0
    size = 40
    case = case
    casefile = casefile
    config_file = CWD + "/" + prefix + "/" + case + "/" + casefile
    odir = "crbdb/$name/sample"
    use_dual_norm = True

    config = g.generateBasisConfig(dim, config_file, time_dependant, odir, case, algo, size, use_dual_norm)

    # compute and save the basis
    g.generate_basis(config=config)



@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_greedy(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp
    dim = dim
    time_dependant = time_dependent
    algo = 1
    size = 60
    case = case
    casefile = casefile
    config_file = CWD + "/" + prefix + "/" + case + "/" + casefile
    odir = "crbdb/$name/greedy"
    tol = 1e-3
    use_dual_norm = True

    config = g.generateBasisConfig(dim, config_file, time_dependant, odir, case, algo, size, tol, use_dual_norm)

    # compute and save the basis
    g.generate_basis(config=config)


@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_POD(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp

    dim = dim
    time_dependant = time_dependent
    algo = 2
    size = 15
    case = case
    casefile = casefile
    config_file = CWD + "/" + prefix + "/" + case + "/" + casefile
    odir = "crbdb/$name/POD"
    tol = 1e-3
    use_dual_norm = True

    config = g.generateBasisConfig(dim, config_file, time_dependant, odir, case, algo, size, tol, use_dual_norm)


    # compute and save the basis
    g.generate_basis(config=config)