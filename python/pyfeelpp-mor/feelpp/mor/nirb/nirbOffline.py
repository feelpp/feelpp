import sys, os
import feelpp
from feelpp.mor.nirb import nirbOffline, nirbOnline, ToolboxModel
from feelpp.mor.nirb.greedy import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
from feelpp.mor.nirb.nirb_perf import checkConvergence
import time
import json
import argparse
import pandas as pd


default_values = {
    "nbSnapshots": 10,
    "outdir": None,

    "doRectification": True,
    "doBiorthonormal": False,
    "doGreedy": True,

    "greedy-Ninit": 5,
    "greedy-Ntrain": 1000,
    "greedy-tol": 1e-5,
    "greedy-Nmax": 50,
    "greedy-samplingMode": "log-random",
    "greedy-Mexpr": "N-1",
}

def get_config_parameter(config, parameter):
    """
    Get parameter from config file, if not found, return default value

    Parameters
    ----------
    config : dict
        configuration of the NIRB model
    parameter : str
        parameter to get from config file

    Returns
    -------
    parameter
    """
    if parameter.startswith("greedy-"):
        return config["greedy-parameters"].get(parameter[7:], default_values[parameter])
    return config.get(parameter, default_values[parameter])


def run_offline_pod(config_nirb):
    """Run offline step using POD

    Parameters
    ----------
    config_nirb : dict
        configuration of the NIRB model

    Returns
    -------
    nirbOffline
        nirbOffline object, with the basis initialized
    """
    nirb_off = nirbOffline(**config_nirb, initCoarse=True)
    nirb_off.initModel()

    nirb_off.generateOperators(coarse=True)
    nbSnap = config_nirb['nbSnapshots']
    _ = nirb_off.computeSnapshots(nbSnap)
    RIC = nirb_off.generateReducedBasis()
    return nirb_off



def run_offline_greedy(config_nirb, Ninit, Ntrain, tol=1e-5, Xi_train=None, Nmax=50, samplingMode="log-random", Mexpr="N-1"):
    """Run offline step using greedy algorithm

    Parameters
    ----------
    config_nirb : dict
        configuration of the NIRB model
    Ninit : int
        size of the initial basis computed using POD
    Ntrain : int
        size of the training set
    tol : float, optional
        tolerance for the greedy algorithm
    Xi_train : list, optional
        training set, by default None. If None, a sample of size Ntrain is generated
    Nmax : int, optional
        maximal size of the basis, by default 50
    samplingMode : str, optional
        mode for sampling the training set, by default "log-random"
    Mexpr : str, optional
        evaluation of M regarding N

    Returns
    -------
    nirbOffline
        nirbOffline object, with the basis initialized
    """

    tb = ToolboxModel(**config_nirb)
    tb.initModel(initCoarse=True)
    interpolator = tb.createInterpolator(tb.tbCoarse, tb.tbFine)

    nirb_off = nirbOffline(initCoarse=True, **config_nirb)
    nirb_off.setModel(tb)

    nirb_on = nirbOnline(**config_nirb)
    nirb_on.setModel(tb, interpolationOperator=interpolator)

    nirb_off.generateOperators(coarse=True)
    res = initProblemGreedy(nirb_off, nirb_on, Ninit, Ntrain, tol=tol,
                Xi_train=Xi_train, Nmax=Nmax, samplingMode=samplingMode, Mexpr=Mexpr)

    return nirb_off, nirb_on, res


if __name__ == "__main__":
    PWD = os.getcwd()

    parser = argparse.ArgumentParser(description='NIRB Offline')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=10)
    parser.add_argument("--convergence", help="Wether get convergence error [default=0]", type=int, default=0)
    parser.add_argument("--load-sampling", help="Load a sampling previously generated [default=0]", type=int, default=0)
    parser.add_argument("--time", help="Wether to solve stationary problem or not [default=0]", type=int, default=0)

    parser.add_argument("--generate-sampling", help="Generate and save a sampling of given size [default=0]\
        If called, no basis is generated", type=int, default=0)
    parser.add_argument("--sampling-path", help="Path to sampling file, for loading or saving [default=\"sampling_train.sample\"]", type=str, default="sampling_train.sample")
    parser.add_argument("--sampling-mode", help="Mode for sampling the training set [default=\"log-random\"]", type=str, default="log-random")
    parser.add_argument("--print-default-values", help="Print default values for config file [default=0]", type=int, default=0)

    args = parser.parse_args()
    config_file = args.config_file

    sampling_path = args.sampling_path

    if args.print_default_values:
        if feelpp.Environment.isMasterRank():
            print("[NIRB] Default values:")
            for k, v in default_values.items():
                print(f"{k:>24}: {v}")
        sys.exit(0)

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']
    outdir = get_config_parameter(config_nirb, 'outdir')

    convergence = args.convergence != 0
    nbSnap = get_config_parameter(config_nirb, 'nbSnapshots')

    doGreedy = get_config_parameter(config_nirb, 'doGreedy')
    doRectification = get_config_parameter(config_nirb, 'doRectification')
    doBiorthonormal = get_config_parameter(config_nirb, 'doBiorthonormal')
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]
    if outdir is None:
        RESPATH = os.path.join(os.getcwd(), "results", rectPath, greedyPath)
    else:
        RESPATH = outdir

    # Create the directory if it does not exist to avoid errors
    if not os.path.exists(RESPATH):
        os.makedirs(RESPATH)

    Xi_train = None
    Xi_train_path = os.path.join(RESPATH, sampling_path)

    ###
    # Only once: generate and save a sampling
    if args.generate_sampling != 0:
        nirb_off = nirbOffline(**config_nirb, initCoarse=False)
        nirb_off.initModel()
        generatedAndSaveSampling(nirb_off.Dmu, args.generate_sampling, path=Xi_train_path, samplingMode=args.sampling_mode)
        sys.exit(0)

    ###
    # If wanted: load the savec sampling to use it in algorithm generation
    if args.load_sampling != 0:
        s = nirb_off.Dmu.sampling()
        N_train = s.readFromFile(Xi_train_path)
        Xi_train = s.getVector()

    # Initialize objects and run the offline stage
    start = time.time()
    if doGreedy:
        if feelpp.Environment.isMasterRank():
            print("+------------------------------------------------+\n[NIRB] Running offline greedy")
        Ninit = get_config_parameter(config_nirb, 'greedy-Ninit')
        Ntrain = get_config_parameter(config_nirb, 'greedy-Ntrain')
        tol = get_config_parameter(config_nirb, 'greedy-tol')
        Nmax = get_config_parameter(config_nirb, 'greedy-Nmax')
        samplingMode = get_config_parameter(config_nirb, 'greedy-samplingMode')
        Mexpr = get_config_parameter(config_nirb, 'greedy-Mexpr')
        nirb_off, nirb_on, _ = run_offline_greedy(config_nirb, Ninit, Ntrain, tol=tol,
                Xi_train=Xi_train, Nmax=Nmax, samplingMode=samplingMode, Mexpr=Mexpr)
    else:
        if feelpp.Environment.isMasterRank():
            print("+------------------------------------------------+\n[NIRB] Running offline")
        nirb_off = run_offline_pod(config_nirb)

    tolortho = 1.e-8

    # nirb_off.orthonormalizeL2(tol=tolortho)
    # nirb_off.orthonormalizeH1(tol=tolortho)
    # nirb_off.orthonormalizeL2(tol=tolortho)

    nirb_off.saveData(RESPATH, force=True)

    finish = time.time()

    print(f"proc {nirb_off.worldcomm.localRank()} Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized(tol=tolortho))
    print(f"proc {nirb_off.worldcomm.localRank()} Is H1 orthonormalized ?", nirb_off.checkH1Orthonormalized())

    if convergence:
        nirb_on = nirb_on = nirbOnline(**config_nirb)
        nirb_on.setModelFromOffline(nirb_off)
        nirb_on.setBasis(nirb_off)
        print(f"[NIRB] Check convergence with sampling [{sampling_path}]")
        Xi_test_path = os.path.join(RESPATH, sampling_path)
        if os.path.isfile(Xi_test_path):
            s = nirb_off.Dmu.sampling()
            N = s.readFromFile(Xi_test_path)
            Xi_test = s.getVector()
            if nirb_off.worldcomm.isMasterRank():
                print(f"[NIRB] Xi_test loaded from {Xi_test_path}")
        else :
            Ns = 30
            Xi_test = generatedAndSaveSampling(nirb_off.Dmu, Ns, path=Xi_test_path, samplingMode="log-random")

        Err = checkConvergence(nirb_off, nirb_on, Ns=30, Xi_test=Xi_test)
        # Err = nirb_off.checkConvergence(Ns=30, Xi_test=Xi_test)
        df = pd.DataFrame(Err)
        errorfile = os.path.join(RESPATH, "offlineError.csv")
        df.to_csv(errorfile, index=False)
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Offline error saved in {os.path.join(os.getcwd(), errorfile)}")


    perf = []
    perf.append(nbSnap)
    perf.append(finish - start)

    if doRectification:
        perfFile = os.path.join(RESPATH, f'nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}_rectif.dat')
    else :
        perfFile = os.path.join(RESPATH, f'nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}.dat')
    WriteVecAppend(perfFile, perf)

    info = nirb_off.getOfflineInfos()

    if nirb_off.worldcomm.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Perf saved in {perfFile}")
        print(f"[NIRB] Offline Elapsed time = ", finish - start)
        print(f"[NIRB] doRectification : {nirb_off.doRectification}, doGreedy : {nirb_off.doGreedy}")
        print(f"[NIRB] Offline part Done !")
