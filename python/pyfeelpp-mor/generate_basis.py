import feelpp.mor.reducedbasis.reducedbasis as mor_rb
import sys, os
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp


import argparse

HOME = os.environ['HOME']

parser = argparse.ArgumentParser()
parser.add_argument('--config-file', help="path to cfg file", type=str)
parser.add_argument('--odir', help="path to output directory", type=str)
parser.add_argument('--dim', help="dimension of the case", type=int)
parser.add_argument('--time-dependant', help="time dependand case", type=bool, default=False)
parser.add_argument('--algo', help="compute using greedy algorithm (default 1) 0 : From sample, 1 : Greedy, 2 : POD", type=int, default=1)
parser.add_argument('--train-size', help="size of the (random) training set", type=int, default=40)
parser.add_argument('--tol', help="tolerance for generating", type=float, default=1e-6)

algos_names = ["generation_from_sample", "greedy", "POD_modes"]

if __name__ == '__main__':
    args = parser.parse_args()
    dim = args.dim
    config_file = args.config_file
    time_dependant = args.time_dependant
    odir = args.odir
    algo = args.algo
    size = args.train_size
    tol = args.tol
else:
    dim = None
    config_file = None
    time_dependant = None
    odir = None
    compute_greedy = None
    algo = None
    size = None
    tol = 1e-6


if config_file is not None:
    split = config_file.split('/')
    case = '/'.join(split[:-1])
    if odir is None:
        dir = HOME + '/feel/reducedbasis/'+split[-2]
    else:
        dir = odir
    casefile = split[-1]


def generate_basis(worldComm=None):

    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()

    if dim not in [2,3]:
        raise ValueError("dim must be 2 or 3")
    if config_file is None:
        raise ValueError("invalid path to config-file")
    if algo not in [0,1,2]:
        raise ValueError(f"Algo {algo} not valid")


    if worldComm.isMasterRank():
        print("==============================================")
        print("Generation of the reduced basis")
        print("           Config-file :", f"{case}/{casefile}")
        print("Data will be stored in :", dir)
        print("Current working directory is ", os.getcwd())
        print("===============================================")


    feelpp.Environment.setConfigFile(f'{case}/{casefile}')
    feelpp.Environment.changeRepository(directory=f'{dir}')

    # Set the toolboxes
    heatBox = heat(dim=dim, order=1)
    heatBox.init()

    modelParameters = heatBox.modelProperties().parameters()
    default_parameter = modelParameters.toParameterValues()


    # model = toolboxmor_2d() if DIM==2 else toolboxmor_3d()
    model = toolboxmor(dim=dim, time_dependent=time_dependant)
    model.setFunctionSpaces( Vh=heatBox.spaceTemperature() )


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

    heatBoxMDEIM = heat(dim=dim, order=1)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    def assembleOnlineMDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)
    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    Dmu = model.parameterSpace()


    def listOfParams(n):
        mus = []
        for _ in range(n):
            mus.append(Dmu.element(True, True))
        return mus

    mubar = Dmu.element(True, False)
    mubar.setParameters(default_parameter)
    if worldComm.isMasterRank():
        print("mubar =")
        mubar.view()

    affineDecomposition = model.getAffineDecomposition()
    Aq = affineDecomposition[0]
    Fq = affineDecomposition[1]

    rb = mor_rb.reducedbasisOffline(mor_rb.convertToPetscMat(Aq[0]),
                                    mor_rb.convertToPetscVec(Fq[0][0]),
                                    model, mubar)
    rb.setVerbose(False)
    if worldComm.isMasterRank():
        print("Size of the big problem :", rb.NN)


    mus = listOfParams(size)
    if worldComm.isMasterRank():
        print("[generate_basis] Start generation of the basis using algo", algos_names[algo])
    if algo == 0:
        rb.computeOfflineReducedBasis(mus)

    elif algo == 1:
        mu_0 = Dmu.element()
        rb.greedy(mu_0, mus, eps_tol=tol)

    elif algo == 2:
        rb.generatePOD(mus, eps_tol=tol)
    else:
        pass

    if worldComm.isMasterRank():
        print("[generate_basis] basis generated ! Now computing errors")

    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()

    if worldComm.isMasterRank():
        print("[generate_basis] Done !")

        rb.saveReducedBasis(dir, force=True)

    return rb



if __name__ == '__main__':
    config = feelpp.globalRepository(f'{dir}')
    sys.argv = ['generate-basis']
    o = toolboxes_options("heat")
    o.add(makeToolboxMorOptions())

    e = feelpp.Environment(sys.argv, opts=o, config=config)

    generate_basis()