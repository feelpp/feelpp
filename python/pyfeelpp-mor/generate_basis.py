import feelpp.mor.reducedbasis.reducedbasis_timeOffline as mor_rb
import sys, os
from feelpp.toolboxes.heat import *
import feelpp.toolboxes.core as core
from feelpp.operators import mass
import feelpp.mor as mor
import feelpp
import numpy as np

import argparse

HOME = os.environ['HOME']

parser = argparse.ArgumentParser()
parser.add_argument('--config-file', help="path to cfg file", type=str)
parser.add_argument('--odir', help="path to output directory", type=str)
parser.add_argument('--case', help="name of the case", type=str)
parser.add_argument('--dim', help="dimension of the case", type=int)
parser.add_argument('--time-dependant', help="time dependand case", type=int, default=0)
parser.add_argument('--algo', help="compute using greedy algorithm (default 1) 0 : From sample, 1 : Greedy, 2 : POD",
    type=int, default=1)
parser.add_argument('--train-size', help="size of the (random) training set", type=int, default=40)
parser.add_argument('--tol', help="tolerance for generating", type=float, default=1e-6)
parser.add_argument('--use-dual-norm', help="use dual norm for the error bound", type=int, default=0)
parser.add_argument('--param', help="parameter to use", type=str)

algos_names = ["generation_from_sample", "greedy", "POD_modes"]

class generateBasisConfig():

    def __init__(self, dim=None, config_file=None, time_dependant=False, odir=None,
                 case=None, algo=1, size=40, tol=1e-6, use_dual_norm=False, param=None):
        self.dim = dim
        self.config_file = config_file
        self.time_dependant = time_dependant
        self.odir = odir
        self.case = case
        self.algo = algo
        self.size = size
        self.tol = tol
        self.use_dual_norm = use_dual_norm

        if param is not None:
            self.param = feelpp.read_json(param)
        else:
            self.param = None

    def getMusk(self, model):
        if self.param is None:
            return None
        
        assert len(self.param['ks']) == len(self.param['mus'])
        n = len(self.param['ks'])
        musk = {}
        for i in range(n):
            mu = model.parameterSpace().element()
            mu.setParameters(self.param['mus'][i])
            musk[mu] = self.param['ks'][i]
        return musk




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
        feelpp.mor.reducedbasis.reducedbasis.reducedbasisOffline: Object dealing with the reduced basis.
        The content is also saved on disk
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

    name = config.case.replace("/","-")# + "-np_" +  str(feelpp.Environment.numberOfProcessors())
    name = name.replace('--', '-')
    if "$name" in config.odir:
        config.odir = config.odir.replace("$name", name)


    if worldComm.isMasterRank():
        cwd = os.getcwd()
        config_file = f"{config.config_file}"
        store = feelpp.Environment.rootRepository()+"/"+config.odir
        box_size = max(33, 26+len(cwd), 26+len(config_file), 30+len(store) ) + 1
        print( "╔"+box_size * "═"+"╗")
        print( "║ Generation of the reduced basis".ljust(box_size) + " ║")
        print(f"║            Config-file : {config_file}".ljust(box_size) + " ║")
        print(f"║ Data will be stored in : {store}".ljust(box_size) + " ║")
        print(f"║ Current working directory is {cwd}".ljust(box_size) + " ║")
        print( "╚"+box_size * "═"+"╝")


    feelpp.Environment.setConfigFile(f'{config.config_file}')

    model_path = "$cfgdir/"+os.path.splitext(os.path.basename(config.config_file))[0] + ".json"
    j = feelpp.read_json(model_path)
    try:
        j.pop('PostProcess')
    except KeyError as e:
        print(f"There was no section {e} in the model")


    crb_model_properties = mor.CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
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

    model = mor.toolboxmor(name=name, dim=config.dim, time_dependent=config.time_dependant)
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
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    heatBoxDEIM = heat(dim=config.dim, order=1)
    heatBoxDEIM.setModelProperties(j)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    heatBoxMDEIM = heat(dim=config.dim, order=1)
    heatBoxMDEIM.setModelProperties(j)
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


    mubar = Dmu.element(True, False)
    mubar.setParameters(default_parameter)
    if worldComm.isMasterRank():
        print("mubar =")
        mubar.view()

    affineDecomposition = model.getAffineDecomposition()
    assert len(affineDecomposition) == [2, 3][config.time_dependant]
    Aq_ = affineDecomposition[0]
    Fq_ = affineDecomposition[1]

    Aq = mor_rb.convertToPetscMat(Aq_[0])
    Fq = []
    for f in Fq_:
        Fq.append(mor_rb.convertToPetscVec(f[0]))

    if not config.time_dependant:

        rb = mor_rb.reducedbasisOffline(Aq=Aq, Fq=Fq, model=model, mubar=mubar,
                output_names=output_names, use_dual_norm=config.use_dual_norm)

        rb.setVerbose(False)
        if worldComm.isMasterRank():
            print("Size of the finite element problem :", rb.NN)

        s = Dmu.sampling()
        s.sampling(config.size, "random")
        mus = s.getVector()
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


    else:
        M = mass(trial=heatBox.spaceTemperature(), test=heatBox.spaceTemperature(), range=feelpp.elements(heatBox.mesh()))
        M = M.to_petsc().mat()
        M.assemble()

        tf = config.param['tf']
        K = config.param['K']
        
        rb = mor_rb.reducedbasisTimeOffline(Aq=Aq, Fq=Fq, Mr=[M], model=model, mubar=mubar,
                output_names=output_names, use_dual_norm=config.use_dual_norm, tf=tf, K=K)

        rb.setVerbose(False)
        if worldComm.isMasterRank():
            print("Size of the finite element problem :", rb.NN)


        musk = config.getMusk(model)
        if worldComm.isMasterRank():
            print("[generate_basis] Start generation of the basis")

        rb.generateBasis(musk)

        if worldComm.isMasterRank():
            print("[generate_basis] basis generated ! Now computing errors")

        def g(k): return 1
        g = np.vectorize(g)

        ts = np.linspace(0, tf, K+1)
        glist = g(ts)

        rb.computeOfflineErrorRhs()
        rb.computeOfflineError(glist)

    if worldComm.isMasterRank():
        print(f"[generate_basis] Done ! Size of the reduced problem : {rb.N}")

    s = rb.saveReducedBasis(feelpp.Environment.rootRepository()+"/"+config.odir, force=True)

    return rb



if __name__ == '__main__':
    # config = feelpp.globalRepository(f'{dir}')
    if 'generate_basis.py' in sys.argv[0]:
        sys.argv = sys.argv[1:]

    args = parser.parse_args(sys.argv)
    dim = args.dim
    config_file = args.config_file
    time_dependant = args.time_dependant == 1
    odir = args.odir if args.odir is not None else "$name"
    algo = args.algo
    size = args.train_size
    tol = args.tol
    split = config_file.split('/')
    case = args.case if args.case is not None else "generate_basis"
    use_dual_norm = args.use_dual_norm == 1
    param = args.param if args.param is not None else None

    config = feelpp.globalRepository(f"generate_basis-{case}")
    sys.argv = [f'generate-basis-{case}']
    o = core.toolboxes_options("heat")
    o.add(mor.makeToolboxMorOptions())

    e = feelpp.Environment(sys.argv, opts=o, config=config)


    config = generateBasisConfig(dim, config_file, time_dependant, odir, case, algo, size, tol, use_dual_norm, param)

    generate_basis(config=config)
