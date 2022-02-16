import shutil, os
import pytest
import pandas as pd

FEEL_PATH = os.environ['HOME'] + '/feel'

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp

from feelpp.mor.reducedbasis.reducedbasis import *

#        (( prefix, case, casefile, dim, use_cache, time_dependant), name     )
cases = [
        #  (('thermal-fin', '2d', 'thermal-fin.cfg', 2, False, False), 'thermal-fin-2d'),
         (('thermal-fin', '2d', 'thermal-fin.cfg', 2, True, False), 'thermal-fin-2d-cached'),
         (('thermal-fin', '3d', 'thermal-fin.cfg', 3, False, False), 'thermal-fin-3d'),
         (('thermal-fin', '3d', 'thermal-fin.cfg', 3, True, False), 'thermal-fin-3d-cached')
        ]
cases_params, cases_ids = list(zip(*cases))



def init_toolbox(prefix, case, casefile, dim, use_cache):

    if not use_cache:
        print('removing cache...')
        try:
            shutil.rmtree(FEEL_PATH + '/crbdb')
            shutil.rmtree(FEEL_PATH + f'{prefix}/{case}')
        except FileNotFoundError:
            print('You asked to remove cache, but there was none, so no problem !')

    feelpp.Environment.setConfigFile(f'{prefix}/{case}/{casefile}')
    feelpp.Environment.changeRepository(directory=f'{prefix}/{case}')

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    return heatBox, dim



def init_model(prefix, case, casefile, dim, use_cache, time_dependent):
    heatBox, dim = init_toolbox(prefix, case, casefile, dim, use_cache)
    model = toolboxmor(dim=dim, time_dependent=time_dependent)

    Dmu = model.parameterSpace()
    mubar = Dmu.element()

    modelProperties = heatBox.modelProperties()
    param = modelProperties.parameters()
    for p in param:
        assert(p[1].hasMinMax())    # check that parameters can vary
        mubar.setParameterNamed(p[1].name(), p[1].value())

    model.setFunctionSpaces(Vh=heatBox.spaceTemperature())

    def assembleDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleRhs()

    def assembleMDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleMatrix()

    model.setAssembleDEIM(fct=assembleDEIM)


    model.setAssembleDEIM(fct=assembleDEIM)
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    heatBoxDEIM = heat(dim=dim, order=1)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0, mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)

    heatBoxMDEIM = heat(dim=dim, order=1)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineMDEIM(mu):
        for i in range(0, mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    return heatBox, model, time_dependent, mubar, assembleMDEIM, assembleDEIM


def init_environment(prefix, case, casefile, dim, use_cache, time_dependent):
    heatBox, model, time_dependent, mubar, assembleMDEIM, assembleDEIM = init_model(prefix, case, casefile, dim, use_cache, time_dependent)

    decomposition = model.getAffineDecomposition()
    assert len(decomposition) == [2,3][time_dependent]

    return heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM


def compar_from_sampling(rb, xi_test):
    """Compare the solutions on the generating sample (should be computer 0)
       xi_test has to be equal to xi_train used to generate the basis
    """
    for mu in xi_test:
        relErrOnU = rb.compareSols(mu)
        assert( relErrOnU < 1e-12 )



def compar_matrix(rb, assembleMDEIM):
    """Compares the construction of the matrix to the toolbox one
    """
    model = rb.model
    Dmu = model.parameterSpace()
    
    mu = Dmu.element(True, False)

    beta = model.computeBetaQm(mu)
    betaA = beta[0]
    M_tb = assembleMDEIM(mu).mat()
    M_tb.assemble()
    M_rb = rb.assembleA(betaA[0])

    assert(M_tb.size == M_rb.size)
    norm = (M_tb - M_rb).norm() / M_tb.norm()
    print(f"relErr = {norm}\n||M_tb|| = {M_tb.norm()}, ||M_rb|| = {M_rb.norm()}")
    assert (norm < 1e-10)


def compar_rhs(rb, assembleDEIM):
    """Compares the construction of the rhs to the toolbox one
    """
    model = rb.model
    Dmu = model.parameterSpace()

    mu = Dmu.element(True, False)
    beta = model.computeBetaQm(mu)
    assert len(beta) == 2, f"len(beta)={len(beta)} and should be 2"
    betaF = beta[1]
    F_tb = assembleDEIM(mu).vec()
    F_tb.assemble()
    F_rb = rb.assembleF(betaF[0][0])
    assert F_tb.size == F_rb.size, f"F_tb = {F_tb.size} != {F_rb.size} = F_rb"
    
    norm = (F_tb-F_rb).norm() / F_tb.norm()
    print(f"relErr = {norm}\n||F_tb|| = {F_tb.norm()}, ||F_rb|| = {F_rb.norm()}")
    assert norm < 1e-10, f"relative error {norm} too high"

def compar_sols(rb, assembleMDEIM, heatBox):
    """Compares the computation of the output
    """
    model = rb.model
    Dmu = model.parameterSpace()
    mu = Dmu.element(True, False)
    # mu.setParameters({"Bi":0.01, "k_0":1, "k_1":0.1, "k_2":0.1, "k_3":0.1, "k_4":0.1})
    beta = model.computeBetaQm(mu)
    assert(len(beta) == 2)
    assembleMDEIM(mu)
    heatBox.solve()
    s_tb = feelpp.mean(range=feelpp.markedfaces(heatBox.mesh(), "Gamma_root"), expr=heatBox.fieldTemperature())[0]

    _,sN = rb.getSolutions(mu)
    
    norm = (sN - s_tb ) / s_tb
    print(f"relErr = {norm}\n||s_tb|| = {s_tb}, ||s_rb|| = {sN}")
    assert norm < 1e-10, f"relative error {norm} is too high"


def compar_solFE(rb, assembleMDEIM, heatBox):
    """Compares the construction of the matrix to the toolbox one
    """
    model = rb.model
    Dmu = model.parameterSpace()
    mu = Dmu.element(True, False)
    beta = model.computeBetaQm(mu)
    assert(len(beta) == 2)
    assembleMDEIM(mu)
    heatBox.solve()
    u_tb = heatBox.fieldTemperature().to_petsc().vec()

    u_fe,_ = rb.getSolutionsFE(mu)

    assert(u_tb.size == u_fe.size)
    
    norm = (u_tb - u_fe).norm() / u_tb.norm()
    print(f"relErr = {norm}\n||u_tb|| = {u_tb}, ||u_fe|| = {u_fe}")
    assert norm < 1e-10, f"relative error {norm} is too high"


@pytest.mark.parametrize("prefix,case,casefile,dim,use_cache,time_dependent", cases_params, ids=cases_ids)
def test_init_reducedbasis(prefix, case, casefile, dim, use_cache, time_dependent, init_feelpp):
    e = init_feelpp    
    heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM = init_environment(prefix, case, casefile, dim, use_cache, time_dependent)

    Aq = decomposition[0]
    Fq = decomposition[1]

    rb = reducedbasis(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar)

    print("\nCompute basis and orthonormalize it")
    def listOfParams(n):
        res = []
        for i in range(n):
            res.append(model.parameterSpace().element())
        return res

    mus = listOfParams(40)
    rb.computeOfflineReducedBasis(mus, orth=True)
    assert( rb.test_orth() )


    print("\nCompute offline error")
    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()

    if size == 1:
        print("\nError with parameters from sampling")
        compar_from_sampling(rb, mus)

    print("\nCompar matrix")
    compar_matrix(rb, assembleMDEIM)

    print("\nCompar Rhs")
    compar_rhs(rb, assembleDEIM)

    print("\nCompar solutions")
    compar_sols(rb, assembleDEIM, heatBox)

    print("\nCompar FE solutions")
    compar_solFE(rb, assembleDEIM, heatBox)








# # Greedy Tests

# @pytest.mark.long
# # @pytest.mark.dependency(depends=['test_init_environment'])
# def test_runGreedy():
#     """runs the greedy algorithm to generate a basis
#        (this test is quite long)
#     """
#     Aq = pytest.decomposition[0]
#     Fq = pytest.decomposition[1]
#     pytest.rbGreedy = reducedbasis(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar, alphaLB)

#     Xi_train = listOfParams(100)
#     mu0 = pytest.Dmu.element(True, True)
#     S = pytest.rbGreedy.greedy(mu0, Xi_train, Nmax=60)

#     assert( pytest.rbGreedy.DeltaMax[-1] < 1e-6 )

# # @pytest.mark.greedy
# # @pytest.mark.dependency(depends=['test_runGreedy'])
# # def test_cvgError(xi_test=None):
# #     if xi_test is None:
# #         xi_test = listOfParams(50)
# #     nb = len(xi_test)

# #     df = pd.DataFrame(columns=['minS', 'maxS', 'meanS', 'minU', 'maxU', 'meanU'],
# #                       index = list(range(1,pytest.rbGreedy.N+1)))

# #     # print("size minS maxS meanS minU maxU meanU")

# #     for size in range(1,pytest.rbGreedy.N+1):
# #         S = np.zeros(nb)
# #         U = np.zeros(nb)

# #         for i in range(nb):
# #             uN, sN = pytest.rbGreedy.getSolutions(xi_test[i], size=size)
# #             u , s  = pytest.rbGreedy.getSolutionsFE(xi_test[i])

# #             u_proj = pytest.rbGreedy.projFE(uN)

# #             S[i] = np.abs(sN - s)/np.abs(s)
# #             U[i] = pytest.rbGreedy.normA(u - u_proj)/pytest.rbGreedy.normA(u)

# #         df.loc[size] = pd.Series({'minS':np.min(S), 'maxS':np.max(S), 'meanS':np.mean(S),
# #                                   'minU':np.min(U), 'maxU':np.max(U), 'meanU':np.mean(U)})
# #         # print(size, np.min(S), np.max(S), np.mean(S), np.min(U), np.max(U), np.mean(U))
# #     return df





# # Other tests

# # @pytest.mark.dependency(depends=['test_computeBasis'])
# def test_save_load():
#     """Tests that the loaded matrices are identical to the ones saved
#     """
#     os.system("rm -rf /tmp/rb")
#     assert( reducedbasis.loadReducedBasis('/tmp/rb', pytest.model) == None )

#     pytest.rb.saveReducedBasis('/tmp/rb')
#     rbLoaded = reducedbasis.loadReducedBasis('/tmp/rb', pytest.model)

#     assert( pytest.rb.Qa == rbLoaded.Qa )
#     assert( pytest.rb.Qf == rbLoaded.Qf )
#     assert( pytest.rb.N == rbLoaded.N )

#     for q in range(rbLoaded.Qa):
#         assert( (pytest.rb.ANq[q] == rbLoaded.ANq[q]).all() ), f"q = {q}"
#     for p in range(rbLoaded.Qf):
#         assert( (pytest.rb.FNp[p] == rbLoaded.FNp[p]).all() ), f"p = {p}"
