import  os
import pytest
import feelpp.mor.generate_basis as g


#        (( prefix, case, casefile, dim, time_dependant), name     )
cases = [
         (('testcase/thermal-fin', '2d', 'thermal-fin.cfg', 2, False), 'thermal-fin-2d'),
         (('testcase/thermal-fin', '3d', 'thermal-fin.cfg', 3, False), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))

OUTDIR = "/tmp/test_reduced_basis"
os.system(f"rm -rf {OUTDIR}")


@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_sample(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp

    os.system(f"rm -rf {OUTDIR}/sample")
    g.dim = dim
    g.time_dependant = time_dependent
    g.compute_greedy = False
    g.size = 40
    g.case = prefix + "/" + case
    g.casefile = casefile
    g.dir = OUTDIR + "/sample"
    g.config_file = ''

    # compute and save the basis
    g.generate_basis()



@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_compute_basis_greedy(prefix, case, casefile, dim, time_dependent, init_feelpp):

    e = init_feelpp
    os.system(f"rm -rf {OUTDIR}/greedy")
    g.dim = dim
    g.time_dependant = time_dependent
    g.compute_greedy = True
    g.size = 60
    g.case = prefix + "/" + case
    g.casefile = casefile
    g.dir = OUTDIR + "/greedy"
    g.config_file = ''
    g.tol = 1e-3

    # compute and save the basis
    g.generate_basis()