import sys
import pytest
import pandas as pd


from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp.core as fppc

from feelpp.mor.reducedbasis.reducedbasis_time import *


# def test_init_argv(configfile):
#     sys.argv += ['--config-file', configfile]
#     assert False, configfile

# if "configfile" in metafunc.fixturenames:
        # metafunc.parametrize("configfile", metafunc.config.getoption("configfile"))

# sys.argv += ['--config-file', '/home/saigre/Documents/rb/thermal-fin/2d/om.cfg']
sys.argv += ['--config-file', '/home/saigre/Documents/rb/thermal-fin/2d/thermal-fin.cfg']



DIM = 2
assert( DIM in [2,3] )


# Set the environment
o=toolboxes_options("heat")
o.add(makeToolboxMorOptions())
e=fppc.Environment(sys.argv,opts=o)


# Set the toolboxes
# TODO: get DIM and time_dependent from cfg file
heatBox = heat(dim=DIM,order=1)
heatBox.init()
# model = toolboxmor_2d() if DIM==2 else toolboxmor_3d()
model = toolboxmor(dim=DIM, time_dependent=True)
model.setFunctionSpaces(Vh=heatBox.spaceTemperature())



# Offline computations

def assembleDEIM(mu):
    for i in range(0, mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleRhs()

def assembleMDEIM(mu):
    for i in range(0, mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleMatrix()

model.setAssembleDEIM(fct=assembleDEIM)
model.setAssembleMDEIM(fct=assembleMDEIM)
model.initModel()

heatBoxDEIM = heat(dim=DIM, order=1)
meshDEIM = model.getDEIMReducedMesh()
heatBoxDEIM.setMesh(meshDEIM)
heatBoxDEIM.init()

def assembleOnlineDEIM(mu):
    for i in range(0, mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
    heatBoxDEIM.updateParameterValues()
    return heatBoxDEIM.assembleRhs()

model.setOnlineAssembleDEIM(assembleOnlineDEIM)

heatBoxMDEIM = heat(dim=DIM, order=1)
meshMDEIM = model.getMDEIMReducedMesh()
heatBoxMDEIM.setMesh(meshMDEIM)
heatBoxMDEIM.init()

def assembleOnlineMDEIM(mu):
    for i in range(0, mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
    heatBoxMDEIM.updateParameterValues()
    return heatBoxMDEIM.assembleMatrix()

model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

model.postInitModel()
model.setInitialized(True)
Dmu = model.parameterSpace()

print("Element :")
Dmu.element().view()
print("\n\n\n")

def listOfParams(n):
    mus = []
    for _ in range(n):
        mus.append(Dmu.element(True, True))
    return mus


mubar = Dmu.element(True, False)    # TODO : see how to get the values from json
mubar.setParameters({"Bi":0.1, "k_0":1, "k_1":1, "k_2":1, "k_3":1, "k_4":1, "time":0})

def alphaLB(mu):
    return min(mu.parameterNamed("k_1"), 1)


# @pytest.mark.dependency(depends=['test_init_argv'])
@pytest.mark.dependency()
def test_init_environment():
    pytest.decomposition = model.getAffineDecomposition()
    assert(len(pytest.decomposition) == 3)




@pytest.mark.dependency(depends=['test_init_environment'])
def test_init_reducedbasis():
    Aq = pytest.decomposition[0]
    Fq = pytest.decomposition[1]
    Mq = pytest.decomposition[2]
    pytest.rb = reducedbasis_time(convertToPetscMat(Aq[0]), convertToPetscVec(Fq[0][0]), model, mubar, alphaLB, convertToPetscMat(Mq[0]), 10, 100)



@pytest.mark.dependency(depends=['test_init_reducedbasis'])
def test_computeBasis():
    """Checks that the reduced basis is well computed and orthonormalized
    """
    pytest.mus = listOfParams(40)
    # pytest.rb.computeOfflineReducedBasis(pytest.mus, orth=True)
    musk = {}
    default_values = {'k0':1, 'k1':1, 'k2':1, 'k3':1, 'k4':1}
    for Bi,ts in zip([0.01, 0.1, 1],[[1,5,10,20,30],[5,10,20],[5,10]]):
        mu = Dmu.element()
        mu.setParameters(default_values)
        mu.setParameters({'Bi':Bi})
        musk[mu] = ts
    pytest.rb.generateBasis(musk)
    # assert( rbTest.test_orth() == np.eye(rbTest.N)).all()
    assert( pytest.rb.test_orth() )



@pytest.mark.dependency(depends=['test_computeBasis'])
def test_computeOfflineError():
    """Tests thaht the compute offline error is well run
    """
    t = np.linspace(0, pytest.rb.tf, pytest.rb.K+1)
    g = 1 - np.cos(t)
    pytest.rb.computeOfflineErrorRhs()
    pytest.rb.computeOfflineError(g)



@pytest.mark.dependency(depends=['test_init_reducedbasis'])
def test_comparMatrix():
    """Compares the construction of the matrix to the toolbox one
    """
    mu = Dmu.element(True, False)    # TODO : see how to get the values from json
    beta = model.computeBetaQm(mu)
    assert(len(beta) == 3)
    betaA = beta[0]
    M_tb = assembleMDEIM(mu).mat()
    M_tb.assemble()
    M_rb = pytest.rb.assembleA(betaA[0])

    assert(M_tb.size == M_rb.size)
    
    norm = (M_tb - M_rb).norm() / M_tb.norm()
    assert norm < 1e-10, f"relative error {norm} is too high"


@pytest.mark.dependency(depends=['test_init_reducedbasis'])
def test_comparRhs():
    """Compares the construction of the rhs to the toolbox one
    """
    mu = Dmu.element(True, False)    # TODO : see how to get the values from json
    mu.setParameters({"Bi":0.01, "k_0":1, "k_1":0.1, "k_2":0.1, "k_3":0.1, "k_4":0.1})
    beta = model.computeBetaQm(mu)
    assert len(beta) == 3, f"len(beta)={len(beta)} and should be 3"
    betaF = beta[1]
    F_tb = assembleDEIM(mu).vec()
    F_tb.assemble()
    F_rb = pytest.rb.assembleF(betaF[0][0])

    assert F_tb.size == F_rb.size, f"F_tb = {F_tb.size} != {F_rb.size} = F_rb"
    norm = (F_tb-F_rb).norm() / F_tb.norm()
    assert norm < 1e-10, f"relative error {norm} too high"


