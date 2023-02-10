import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.greedy import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import json
import argparse
import pandas as pd



def run_offline(config_nirb, regulParam=1e-10):
    """Run offline step using POD

    Parameters
    ----------
    config_nirb : dict
        configuration of the NIRB model
    regulParam : float, optional
        Regularization parameter, by default 1e-10

    Returns
    -------
    nirbOffline
        nirbOffline object, with the basis initialized
    """
    nirb_off = nirbOffline(**config_nirb, initCoarse=True)
    nirb_off.initModel()

    nirb_off.generateOperators(coarse=True)
    nbSnap = config_nirb['nbSnapshots']
    _ = nirb_off.initProblem(nbSnap)
    RIC = nirb_off.generateReducedBasis(regulParam=regulParam)
    print(f"[NIRB] Basis of size {nirb_off.N} generated, RIC = {RIC}")

    return nirb_off



def run_offline_greedy(config_nirb, Ninit, Ntrain, eps=1e-5, Xi_train=None, Nmax=50, samplingMode="log-random", Mexpr="N-1"):
    """Run offline step using greedy algorithm

    Parameters
    ----------
    config_nirb : dict
        configuration of the NIRB model
    Ninit : int
        size of the initial basis computed using POD
    Ntrain : int
        size of the training set
    eps : float, optional
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
    res = initProblemGreedy(nirb_off, nirb_on, Ninit, Ntrain, eps=eps,
                Xi_train=Xi_train, Nmax=Nmax, samplingMode=samplingMode, Mexpr=Mexpr)

    return nirb_off, nirb_on, res


if __name__ == "__main__":
    PWD = os.getcwd()

    parser = argparse.ArgumentParser(description='NIRB Offline')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=10)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)
    parser.add_argument("--greedy", help="Wether with or without Greedy [default=0]", type=int, default=0)
    parser.add_argument("--biortho", help="Wether with or without bi-orthonormalization [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Wether get convergence error [default=0]", type=int, default=0)

    parser.add_argument("--generate-sampling", help="Generate and save a sampling of given size [default=0]\
        If called, no basis is generated", type=int, default=0)
    parser.add_argument("--load-sampling", help="Load a sampling previously generated [default=0]", type=int, default=0)
    parser.add_argument("--sampling-path", help="Path to sampling file, for loading or saving [default=\"sampling_train.sample\"]", type=str, default="sampling_train.sample")

    # Greedy arguments
    parser.add_argument("--Ninit", help="[Greedy] Number of initial snapshots [default=5]", type=int, default=5)
    parser.add_argument("--Ntrain", help="[Greedy] Number of training snapshots [default=1000]", type=int, default=1000)
    parser.add_argument("--eps", help="[Greedy] Tolerance [default=1e-5]", type=float, default=1e-5)
    parser.add_argument("--Nmax", help="[Greedy] Maximum number of snapshots [default=50]", type=int, default=50)
    parser.add_argument("--sampling-mode", help="[Greedy] Sampling mode [default=log-random]", type=str, default="log-random")
    parser.add_argument("--Mexpr", help="[Greedy] Expression of N as a function of N [default=\"N-1\"]", type=str, default="N-1")

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']

    convergence = args.convergence != 0
    nbSnap = args.N

    config_nirb['greedy-generation'] = args.greedy != 0
    config_nirb['doBiorthonormal'] = args.biortho != 0

    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]

    if outdir is None:
        RESPATH = f"results/{rectPath}/{greedyPath}"
    else:
        RESPATH = outdir

    ###
    # Only once: generate and save a sampling
    if args.generate_sampling != 0:
        nirb_off = nirbOffline(**config_nirb, initCoarse=False)
        nirb_off.initModel()
        generatedAndSaveSampling(nirb_off.Dmu, args.generate_sampling, path=Xi_train_path)
        sys.exit(0)

    ###
    # If wanted: load the savec sampling to use it in algorithm generation
    Xi_train = None
    Xi_train_path = os.path.join(RESPATH, args.sampling_path)
    if args.load_sampling != 0:
        s = nirb_off.Dmu.sampling()
        N_train = s.readFromFile(Xi_train_path)
        Xi_train = s.getVector()

    # Initialize objects and run the offline stage
    start = time.time()
    if config_nirb['greedy-generation']:
        nirb_off, nirb_on, _ = run_offline_greedy(config_nirb, args.Ninit, args.Ntrain, eps=args.eps,
                Xi_train=Xi_train, Nmax=args.Nmax, samplingMode=args.sampling_mode, Mexpr=args.Mexpr)
    else:
        nirb_off = run_offline(config_nirb)

    tolortho = 1.e-8

    # nirb_off.orthonormalizeL2(tol=tolortho)
    # nirb_off.orthonormalizeH1(tol=tolortho)
    # nirb_off.orthonormalizeL2(tol=tolortho)

    nirb_off.saveData(RESPATH, force=True)

    finish = time.time()

    print(f"proc {nirb_off.worldcomm.localRank()} Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized(tol=tolortho))
    print(f"proc {nirb_off.worldcomm.localRank()} Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    if convergence:
        Xi_test_path = os.path.join(RESPATH, "sampling_test.sample")
        if os.path.isfile(Xi_test_path):
            s = nirb_off.Dmu.sampling()
            N = s.readFromFile(Xi_test_path)
            Xi_test = s.getVector()
            if nirb_off.worldcomm.isMasterRank():
                print(f"[NIRB] Xi_test loaded from {Xi_test_path}")
        else :
            Ns = 30
            Xi_test = generatedAndSaveSampling(nirb_off.Dmu, Ns, path=Xi_test_path, samplingMode="log-random")

        Err = nirb_off.checkConvergence(Ns=30, Xi_test=Xi_test)
        df = pd.DataFrame(Err)
        file = os.path.join(RESPATH, "offlineError.csv")
        df.to_csv(file, index=False)
        print(f"[NIRB] Offline error saved in {os.path.join(os.getcwd(), file)}")


    perf = []
    perf.append(nbSnap)
    perf.append(finish-start)

    if doRectification:
        file = os.path.join(RESPATH, f'nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}_rectif.dat')
    else :
        file = os.path.join(RESPATH, f'nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}.dat')
    WriteVecAppend(file, perf)

    info = nirb_off.getOfflineInfos()

    if nirb_off.worldcomm.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] doRectification : {nirb_off.doRectification}, doGreedy : {nirb_off.doGreedy}")
        print(f"[NIRB] Offline part Done !")
