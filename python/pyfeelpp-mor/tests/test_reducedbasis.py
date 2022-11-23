import shutil, os
import pytest
import time


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp

# when reducedOffline is loaded, reducedbasis is automatically loaded
from feelpp.mor.reducedbasis.reducedbasisOffline import *

#        (( prefix, case, casefile, dim, time_dependant), name     )
cases = [
         (('testcase', 'square/2d', 'testcase2d.cfg', 2, False), 'testcase-2d'),
         (('testcase', 'square/3d', 'testcase3d.cfg', 3, False), 'testcase-3d'),
         (('testcase', 'thermal-fin/2d', 'thermal-fin.cfg', 2, False), 'thermal-fin-2d'),
         (('testcase', 'thermal-fin/2d', 'thermal-fin.cfg', 2, True), 'thermal-fin-2d-time'),
         (('testcase', 'thermal-fin/3d', 'thermal-fin.cfg', 3, False), 'thermal-fin-3d'),
        ]
cases_params, cases_ids = list(zip(*cases))



def init_toolbox(prefix, case, casefile, dim):


    feelpp.Environment.setConfigFile(f'{prefix}/{case}/{casefile}')

    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    return heatBox, dim



def init_model(prefix, case, casefile, dim, time_dependent):
    heatBox, dim = init_toolbox(prefix, case, casefile, dim)
    name = case.replace("/","-") + "-np_" +  str(feelpp.Environment.numberOfProcessors())
    model = toolboxmor(name=name, dim=dim, time_dependent=time_dependent)


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
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    Dmu = model.parameterSpace()
    mubar = Dmu.element()

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"
    model_properties = feelpp.ModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    model_properties.setup(model_path)
    model_parameter = model_properties.parameters()
    default_parameter = model_parameter.toParameterValues()

    mubar.setParameters(default_parameter)

    j = feelpp.read_json(model_path.replace("$cfgdir", feelpp.Environment.expand("$cfgdir")))
    try:
        j.pop('PostProcess')
    except KeyError as e:
        print(f"There was no section {e} in the model")
    
    crb_model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    crb_model_properties.setup(model_path)
    crb_model_outputs = crb_model_properties.outputs()

    output_names = []
    for n, _ in crb_model_outputs:
        output_names.append(n)

    heatBoxDEIM = heat(dim=dim, order=1)
    heatBoxDEIM.setModelProperties(j)
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
    heatBoxMDEIM.setModelProperties(j)
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

    return heatBox, model, time_dependent, mubar, assembleMDEIM, assembleDEIM, output_names


def init_environment(prefix, case, casefile, dim, time_dependent):
    heatBox, model, time_dependent, mubar, assembleMDEIM, assembleDEIM, output_names = init_model(prefix, case, casefile, dim, time_dependent)

    decomposition = model.getAffineDecomposition()
    assert len(decomposition) == [2,3][time_dependent]

    return heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM, output_names


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

    heatBox.exportResults()
    meas = heatBox.postProcessMeasures().values()

    meas_names = list(meas.keys())

    _,sN = rb.getSolutionsFE(mu)

    print(rb.N_output)


    for k in range(rb.N_output):

        output = rb.getOutputName(k)
        i = 0
        while output not in meas_names[i]: i += 1
        _, s_eim = rb.getSolutionsFE(mu, k=k)
        s_tb = meas[meas_names[i]]

        norm = abs(s_eim - s_tb ) / abs(s_tb)
        print(f"k={k} relErr = {norm}\n||s_tb|| = {s_tb}, ||s_eim|| = {s_eim}")
        # assert norm < 1e-10, f"relative error {norm} is too high for output {output}" 



def compar_solFE(rb, assembleMDEIM, heatBox):
    """Compares the solution of the FE to the one of the RB
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
    print(f"relErr = {norm}\n||u_tb|| = {u_tb.norm()}, ||u_fe|| = {u_fe.norm()}")
    assert norm < 1e-10, f"relative error {norm} is too high"


def effectivity(rb):
    for _ in range(20):
        mu = rb.model.parameterSpace().element()
        eta = rb.computeEffectivity(mu)
        UB_LB = rb.UB_LB(mu)

        print(1, eta, UB_LB)

        # assert 1 <= eta and eta <= UB_LB


def save_and_load(rb):
    """Tests that the loaded matrices are identical to the ones saved
    """
    os.system("rm -rf /tmp/test_reducedbasis")

    path = rb.saveReducedBasis("/tmp/test_reducedbasis", force=True)

    rbLoaded = reducedbasis(None)
    rbLoaded.loadReducedBasis(path, rb.model)

    for n in rb.mubar.parameterNames():
        assert rb.mubar.parameterNamed(n) == rbLoaded.mubar.parameterNamed(n)

    assert( rb.Qa == rbLoaded.Qa )
    assert( rb.Qf == rbLoaded.Qf )
    assert( rb.N == rbLoaded.N )

    for q in range(rbLoaded.Qa):
        assert( (rb.ANq[q] == rbLoaded.ANq[q]).all() ), f"q = {q}"
    for p in range(rbLoaded.Qf):
        assert( (rb.FNp[p] == rbLoaded.FNp[p]).all() ), f"p = {p}"

    assert( (rb.SS == rbLoaded.SS).all() ), "SS"
    assert( (rb.SL == rbLoaded.SL).all() ), "SL"
    assert( (rb.LL == rbLoaded.LL).all() ), "LL"
    
    if rb.DeltaMax is None:
        assert rbLoaded.DeltaMax is None, "DeltaMax"
    else:
        assert( (rb.DeltaMax == rbLoaded.DeltaMax).all() ), "DeltaMax"

    assert rb.alphaMubar == rbLoaded.alphaMubar, "alphaMubar"
    assert rb.gammaMubar == rbLoaded.gammaMubar, "gammaMubar"




@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_reducedbasis_sample(prefix, case, casefile, dim, time_dependent, init_feelpp):
    e = init_feelpp    
    heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM, output_names = \
        init_environment(prefix, case, casefile, dim, time_dependent)

    Aq = decomposition[0]
    Fq_ = decomposition[1]

    Fq = []
    for f in Fq_:
        Fq.append(convertToPetscVec(f[0]))

    rb = reducedbasisOffline(Aq=convertToPetscMat(Aq[0]), Fq=Fq, model=model, mubar=mubar, output_names=output_names)

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

    # compute_time_for_offline_solution(rb)
    # compute_time_for_online_error_computation(rb)

    if feelpp.Environment.isMasterRank():
        print("\nCompar matrix")
    compar_matrix(rb, assembleMDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar Rhs")
    compar_rhs(rb, assembleDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar solutions")
    compar_sols(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\nCompar FE solutions")
    compar_solFE(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\n Save and reload basis")
    save_and_load(rb)

    if feelpp.Environment.isMasterRank():
        print("\nTest on effectivity")
    effectivity(rb)




@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_reducedbasis_greedy(prefix, case, casefile, dim, time_dependent, init_feelpp):
    e = init_feelpp
    heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM, output_names = \
        init_environment(prefix, case, casefile, dim, time_dependent)

    Aq = decomposition[0]
    Fq_ = decomposition[1]

    Fq = []
    for f in Fq_:
        Fq.append(convertToPetscVec(f[0]))

    rb = reducedbasisOffline(Aq=convertToPetscMat(Aq[0]), Fq=Fq, model=model, mubar=mubar, use_dual_norm=True)

    print("\nCompute basis and orthonormalize it")
    def listOfParams(n):
        res = []
        for i in range(n):
            res.append( model.parameterSpace().element())
        return res

    Xi_train = listOfParams(500)
    mu0 = model.parameterSpace().element()
    rb.greedy(mu0, Xi_train, eps_tol=1e-4)

    assert( rb.test_orth() )
    assert( rb.DeltaMax[-1] < 1e-4 )

    if feelpp.Environment.isMasterRank():
        print("\nCompute offline error")
    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()

    if feelpp.Environment.isMasterRank():
        print("\nCompar matrix")
    compar_matrix(rb, assembleMDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar Rhs")
    compar_rhs(rb, assembleDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar solutions")
    compar_sols(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\nCompar FE solutions")
    compar_solFE(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\n Save and reload basis")
    save_and_load(rb)




@pytest.mark.parametrize("prefix,case,casefile,dim,time_dependent", cases_params, ids=cases_ids)
def test_reducedbasis_pod(prefix, case, casefile, dim, time_dependent, init_feelpp):
    e = init_feelpp    
    heatBox, model, decomposition, mubar, assembleMDEIM, assembleDEIM, output_names = \
        init_environment(prefix, case, casefile, dim, time_dependent)

    Aq = decomposition[0]
    Fq_ = decomposition[1]

    Fq = []
    for f in Fq_:
        Fq.append(convertToPetscVec(f[0]))

    rb = reducedbasisOffline(Aq=convertToPetscMat(Aq[0]), Fq=Fq, model=model, mubar=mubar)

    print("\nCompute basis and orthonormalize it")
    def listOfParams(n):
        res = []
        for i in range(n):
            res.append( model.parameterSpace().element())
        return res

    Xi_train = listOfParams(15)
    rb.generatePOD(Xi_train, eps_tol=1e-3)

    assert( rb.test_orth() )
    assert( rb.DeltaMax[-1] < 1e-3 )


    if feelpp.Environment.isMasterRank():
        print("\nCompute offline error")
    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()


    if feelpp.Environment.isMasterRank():
        print("\nCompar matrix")
    compar_matrix(rb, assembleMDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar Rhs")
    compar_rhs(rb, assembleDEIM)

    if feelpp.Environment.isMasterRank():
        print("\nCompar solutions")
    compar_sols(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\nCompar FE solutions")
    compar_solFE(rb, assembleDEIM, heatBox)

    if feelpp.Environment.isMasterRank():
        print("\n Save and reload basis")
    save_and_load(rb)


def compute_time_for_offline_solution(rb):
    #function rb.getSolutions

    print(f"File is saved in {os.getcwd()}")
    f = open("compute_time_for_offline_solution.dat", "w")
    f.write("function time(s)\n")
    t_computeBetaQm = 0
    t_assembleAN = 0
    t_assembleFN = 0
    t_solve = 0

    N = 1000
    for i in range(N):
        mu = rb.model.parameterSpace().element()
        t0 = time.process_time()
        beta = rb.model.computeBetaQm(mu)
        t1 = time.process_time()
        A_mu = rb.assembleAN(beta[0][0], size=size)
        t2 = time.process_time()
        F_mu = rb.assembleFN(beta[1][0][0], size=size)
        t3 = time.process_time()

        sol = np.linalg.solve(A_mu, F_mu)
        t4 = time.process_time()

        t_computeBetaQm += t1 - t0
        t_assembleAN += t2 - t1
        t_assembleFN += t3 - t2
        t_solve += t4 - t3

    f.write(f"computeBetaQm {t_computeBetaQm/N}\nassembleAN {t_assembleAN/N}\nassembleFN {t_assembleFN/N}\nsolve {t_solve/N}")
    f.close()


def compute_time_for_online_error_computation(rb):
    # function rb.computeOnlineError
    print(f"File is saved in {os.getcwd()}")
    f = open("compute_time_for_online_error_computation.dat", "w")
    f.write("function time(s)\n")
    t_computeBetaQm = 0
    t_assemble = 0
    t_solve = 0
    t_compute_numpy = 0
    t_compute_loop = 0

    N = 1000
    for i in range(N):
        mu = rb.model.parameterSpace().element()

        t0 = time.process_time()
        beta = rb.model.computeBetaQm(mu)
        t1 = time.process_time()
        betaA = beta[0][0]
        betaF = beta[1][0][0]
        A_mu = rb.assembleAN(betaA)
        F_mu = rb.assembleFN(betaF)
        t2 = time.process_time()
        uN = np.linalg.solve(A_mu, F_mu)


        t3 = time.process_time()
        s1_ = betaF @ rb.SS @ betaF    # beta_p*beta_p'*(Sp,Sp')
        s2_ = np.einsum('q,p,n,qpn', betaA, betaF, uN, rb.SL)
        s3_ = np.einsum('q,r,n,m,qnrm', betaA, betaA, uN, uN, rb.LL, optimize=True)
        t4 = time.process_time()

        s1 = 0
        for p in range(rb.Qf):
            for p_ in range(rb.Qf):
                s1 += betaF[p] * betaF[p_] * rb.SS[p,p_]
        s2 = 0
        for p in range(rb.Qf):
            for q in range(rb.Qa):
                for n in range(rb.N):
                    s2 += betaF[p] * betaA[q] * uN[n] * rb.SL[q,p,n]
        s3 = 0
        for q in range(rb.Qa):
            for q_ in range(rb.Qa):
                for n in range(rb.N):
                    for n_ in range(rb.N):
                        s3 += betaA[q] * betaA[q_] * uN[n] * uN[n_] * rb.LL[q,n,q_,n_]
        t5 = time.process_time()
        
        assert((s1_-s1)/s1 < 1e-10)
        assert((s2_-s2)/s2 < 1e-10)
        assert((s3_-s3)/s3 < 1e-10)

        t_computeBetaQm += t1 - t0
        t_assemble += t2 - t1
        t_solve += t3 - t2
        t_compute_numpy += t4 - t3
        t_compute_loop += t5 - t4

    f.write(f"computeBetaQm {t_computeBetaQm/N:e}\nassemble {t_assemble/N:e}\nsolve {t_solve/N:e}\ncomputeNumpy {t_compute_numpy/N:e}\ncomputeLoop {t_compute_loop/N:e}")
    f.close()
