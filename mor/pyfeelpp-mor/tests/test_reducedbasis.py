import sys
import pytest
import pandas as pd


from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp

from reducedbasis.reducedbasis import *

sys.argv += ['--config-file', '/home/saigre/feel/thermal-fin/pb4/thermal-fin.cfg']

DIM = 2
assert( DIM in [2,3] )


# Set the environment
o=toolboxes_options("heat")
o.add(makeToolboxMorOptions())
e=feelpp.Environment(sys.argv,opts=o)


# Set the toolboxes
# TODO: get DIM and time_dependent from cfg file
heatBox=heat(dim=DIM,order=1)
heatBox.init()
# model = toolboxmor_2d() if DIM==2 else toolboxmor_3d()
model = toolboxmor(dim=DIM, time_dependent=False)
model.setFunctionSpaces( Vh=heatBox.spaceTemperature())



# Offline computations

def assembleDEIM(mu):
    for i in range(0,mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleRhs()

def assembleMDEIM(mu):
    for i in range(0,mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleMatrix()

model.setAssembleDEIM(fct=assembleDEIM)
model.setAssembleMDEIM(fct=assembleMDEIM)
model.initModel()

heatBoxDEIM=heat(dim=DIM,order=1)
meshDEIM = model.getDEIMReducedMesh()
heatBoxDEIM.setMesh(meshDEIM)
heatBoxDEIM.init()

def assembleOnlineDEIM(mu):
    for i in range(0,mu.size()):
        heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBoxDEIM.updateParameterValues()
    return heatBoxDEIM.assembleRhs()

model.setOnlineAssembleDEIM(assembleOnlineDEIM)

heatBoxMDEIM=heat(dim=DIM,order=1)
meshMDEIM = model.getMDEIMReducedMesh()
heatBoxMDEIM.setMesh(meshMDEIM)
heatBoxMDEIM.init()

def assembleOnlineMDEIM(mu):
    for i in range(0,mu.size()):
        heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBoxMDEIM.updateParameterValues()
    return heatBoxMDEIM.assembleMatrix()

model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

model.postInitModel()
model.setInitialized(True)
Dmu = model.parameterSpace()


def listOfParams(n):
    mus = []
    for _ in range(n):
        mus.append(Dmu.element(True, True))
    return mus


mubar = Dmu.element(True, False)    # TODO : see how to get the values from json
mubar.setParameters({"Bi":0.1, "k_0":1, "k_1":1, "k_2":1, "k_3":1, "k_4":1})

def alphaLB(mu):
    return min(mu.parameterNamed("k_1"), 1)


@pytest.mark.dependency()
def test_init_environment():
    pytest.decomposition = model.getAffineDecomposition()
    assert(len(pytest.decomposition) == 2)


@pytest.mark.dependency(depends=['test_init_environment'])
def test_init_reducedbasis():

    Aq = pytest.decomposition[0]
    Fq = pytest.decomposition[1]

    pytest.rb = reducedbasis(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar, alphaLB)



@pytest.mark.dependency(depends=['test_init_reducedbasis'])
def test_computeBasis():
    """Cheks that the reduced basis is well computed and orthonormalized
    """
    pytest.mus = listOfParams(40)
    pytest.rb.computeOfflineReducedBasis(pytest.mus, orth=True)
    # assert( rbTest.test_orth() == np.eye(rbTest.N)).all()
    assert( pytest.rb.test_orth() )



@pytest.mark.dependency(depends=['test_computeBasis'])
def test_computeOfflineError():
    """Tests thaht the compute offline error is well run
    """
    pytest.rb.computeOfflineErrorRhs()
    pytest.rb.computeOfflineError()


@pytest.mark.dependency(depends=['test_computeBasis'])
def test_compar():
    """Compare the solutions on the generating sample (should be computer 0)
    """
    for mu in pytest.mus:
        relErrOnU = pytest.rb.compareSols(mu)
        assert( relErrOnU < 1e-12 )



def test_for_param():
    """check that the error on reduced basis for the generating sample is null, when the basis is not orthonormalized
    """
    Aq = pytest.decomposition[0]
    Fq = pytest.decomposition[1]
    
    rbParam = reducedbasis(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar, alphaLB)
    # reduced basis only when RB not orthonormilized
    rbParam.computeOfflineReducedBasis(pytest.mus, orth=False)
    for i,mu in enumerate(pytest.mus):
        print('check RB {} with mu:{}'.format(i,mu))
        [betaA, betaF] = model.computeBetaQm(mu)
        A = rbParam.assembleA(betaA[0])
        F = rbParam.assembleF(betaF[0][0])

        u,_ = pytest.rb.getSolutionsFE(mu)
        uN = np.array(np.zeros((rbParam.N)))
        uN[i]=1
        AN = rbParam.assembleAN(betaA[0])
        FN = rbParam.assembleFN(betaF[0][0])
        # print('FN:',FN)
        # print("_ NN", "N")
        # print("F", u.dot(F), uN.T @ FN)
        # print("A", u.dot(A * u), uN.T @ AN @ uN)
        assert(abs(u.dot(F) - uN.T @ FN)/abs(u.dot(F)) < 1e-10), "abs(u.dot(F) - uN.T @ FN)/abs(u.dot(F)) = {}".format(abs(u.dot(F) - uN.T @ FN)/abs(u.dot(F)))
        assert(abs(u.dot(A * u) - uN.T @ AN @ uN)/abs(u.dot(A * u)) < 1e-10), "abs(u.dot(A * u) - uN.T @ AN @ uN)/abs(u.dot(A * u)) = {}".format(abs(u.dot(A * u) - uN.T @ AN @ uN)/abs(u.dot(A * u)))







# Greedy Tests

@pytest.mark.dependency(depends=['test_init_environment'])
def test_runGreedy():
    Aq = pytest.decomposition[0]
    Fq = pytest.decomposition[1]
    pytest.rbGreedy = reducedbasis(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar, alphaLB)

    Xi_train = listOfParams(100)
    mu0 = Dmu.element(True, True)
    S = pytest.rbGreedy.greedy(mu0, Xi_train, Nmax=60)

    assert( pytest.rbGreedy.DeltaMax[-1] < 1e-6 )



# Other tests

@pytest.mark.dependency(depends=['test_computeBasis'])
def test_save_load():
    """Tests that the loaded matrices are identical to the ones saved
    """
    os.system("rm -rf /tmp/rb")
    assert( reducedbasis.loadReducedBasis('/tmp/rb', model) == None )

    pytest.rb.saveReducedBasis('/tmp/rb')
    rbLoaded = reducedbasis.loadReducedBasis('/tmp/rb', model)

    assert( pytest.rb.Qa == rbLoaded.Qa )
    assert( pytest.rb.Qf == rbLoaded.Qf )
    assert( pytest.rb.N == rbLoaded.N )

    for q in range(rbLoaded.Qa):
        assert( (pytest.rb.ANq[q] == rbLoaded.ANq[q]).all() ), f"q = {q}"
    for p in range(rbLoaded.Qf):
        assert( (pytest.rb.FNp[p] == rbLoaded.FNp[p]).all() ), f"p = {p}"