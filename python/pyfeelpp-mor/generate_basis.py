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
parser.add_argument('--case', help="name of the case", type=str)
parser.add_argument('--dim', help="dimension of the case", type=int)
parser.add_argument('--time-dependant', help="time dependand case", type=bool, default=False)
parser.add_argument('--algo', help="compute using greedy algorithm (default 1) 0 : From sample, 1 : Greedy, 2 : POD", type=int, default=1)
parser.add_argument('--train-size', help="size of the (random) training set", type=int, default=40)
parser.add_argument('--tol', help="tolerance for generating", type=float, default=1e-6)

algos_names = ["generation_from_sample", "greedy", "POD_modes"]

class generateBasisConfig():

    def __init__(self, dim=None, config_file=None, time_dependant=False, odir=None, case=None, algo=1, size=40, tol=1e-6):
        self.dim = dim
        self.config_file = config_file
        self.time_dependant = time_dependant
        self.odir = odir
        self.case = case
        self.algo = algo
        self.size = size
        self.tol = tol



def generate_basis(worldComm=None, config=None):
    """Generates the reduced basis

    Args:
        worldComm (feelpp._core.WorldComm, optional): Pointer to the MPI communicator. Defaults to None.
        config (generateBasisConfig, optional): Configuration for the generation of the basis. Defaults to None.

    Raises:
        ValueError: if no configuration is given
        ValueError: if dimension if not 2 or 3
        ValueError: if no config_file is given
        ValueError: if value for algorithm is not valid (see documentation)

    Returns:
        feelpp.mor.reducedbasis.reducedbasis.reducedbasisOffline: Object dealing with the reduced basis. The content is also saved on disk
    """

    if config is None:
        raise ValueError("No configuration passed. Please give an instance of generateBasisConfig")

    if worldComm is None:
        worldComm = feelpp.Environment.worldCommPtr()

    if config.dim not in [2,3]:
        raise ValueError("dim must be 2 or 3")
    if config.config_file is None:
        raise ValueError("invalid path to config-file")
    if config.algo not in [0,1,2]:
        raise ValueError(f"Algo {config.algo} not valid")

    name = config.case.replace("/","-") + "-" + "-np_" +  str(feelpp.Environment.numberOfProcessors())
    name = name.replace('--', '-')
    if "$name" in config.odir:
        config.odir = config.odir.replace("$name", name)


    if worldComm.isMasterRank():
        print("==============================================")
        print("Generation of the reduced basis")
        print("           Config-file :", f"{config.config_file}")
        print("Data will be stored in :", feelpp.Environment.rootRepository()+"/"+config.odir)
        print("Current working directory is ", os.getcwd())
        print("===============================================")


    feelpp.Environment.setConfigFile(f'{config.config_file}')

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(config.config_file))[0] + ".json"

    crb_model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
    crb_model_properties.setup(model_path)
    crb_model_outputs = crb_model_properties.outputs()

    output_names = []
    for n, _ in crb_model_outputs:
        output_names.append(n)
    
    if worldComm.isMasterRank():
        print(f"[generatebasis] Outputs of the models are {output_names}")


    # Set the toolboxes
    heatBox = heat(dim=config.dim, order=1)
    heatBox.init()

    modelParameters = heatBox.modelProperties().parameters()
    default_parameter = modelParameters.toParameterValues()

    model = toolboxmor(name=name, dim=config.dim, time_dependent=config.time_dependant)
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

    heatBoxDEIM = heat(dim=config.dim, order=1)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    heatBoxMDEIM = heat(dim=config.dim, order=1)
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
    Fq_ = affineDecomposition[1]

    Fq = []
    for f in Fq_:
        Fq.append(mor_rb.convertToPetscVec(f[0]))

    rb = mor_rb.reducedbasisOffline(mor_rb.convertToPetscMat(Aq[0]), Fq, model, mubar, output_names=output_names)
    rb.setVerbose(False)
    if worldComm.isMasterRank():
        print("Size of the big problem :", rb.NN)


    mus = listOfParams(config.size)
    if worldComm.isMasterRank():
        print("[generate_basis] Start generation of the basis using algo", algos_names[config.algo])
    if config.algo == 0:
        rb.computeOfflineReducedBasis(mus)

    elif config.algo == 1:
        mu_0 = Dmu.element()
        rb.greedy(mu_0, mus, eps_tol=config.tol)

    elif config.algo == 2:
        rb.generatePOD(mus, eps_tol=config.tol)
    else:
        pass

    if worldComm.isMasterRank():
        print("[generate_basis] basis generated ! Now computing errors")

    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()

    if worldComm.isMasterRank():
        print("[generate_basis] Done !")

    s = rb.saveReducedBasis(feelpp.Environment.rootRepository()+"/"+config.odir, force=True)

    return rb



if __name__ == '__main__':
    # config = feelpp.globalRepository(f'{dir}')
    if 'generate_basis.py' in sys.argv[0]:
        sys.argv = sys.argv[1:]
    args = parser.parse_args(sys.argv)
    sys.argv = ['generate-basis']
    o = toolboxes_options("heat")
    o.add(makeToolboxMorOptions())

    e = feelpp.Environment(sys.argv, opts=o)

    dim = args.dim
    config_file = args.config_file
    time_dependant = args.time_dependant
    odir = args.odir if args.odir is not None else "$name"
    algo = args.algo
    size = args.train_size
    tol = args.tol
    split = config_file.split('/')
    case = args.case if args.case is not None else "generate_basis"

    config = generateBasisConfig(dim, config_file, time_dependant, odir, case, algo, size, tol)

    generate_basis(config=config)