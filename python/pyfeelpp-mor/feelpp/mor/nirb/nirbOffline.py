import sys
import feelpp.core as fppc
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import json
import argparse
import pandas as pd


if __name__ == "__main__":
    PWD = os.getcwd()

    parser = argparse.ArgumentParser(description='NIRB Offline')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=10)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)
    parser.add_argument("--greedy", help="Wether with or without Greedy [default=0]", type=int, default=0)
    parser.add_argument("--biortho", help="Wether with or without bi-orthonormalization [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Wether get convergence error [default=0]", type=int, default=0)

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = fppc.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = fppc.Environment.expand(cfg['nirb']['filename'])
    config_nirb = fppc.readJson(nirb_file)['nirb']

    greedy = args.greedy
    conv = args.convergence
    nbSnap=args.N
    doBiortho = args.biortho

    bo = [False, True]
    convergence = bo[conv]


    config_nirb['greedy-generation'] = bo[greedy]
    config_nirb['doBiorthonormal'] = bo[doBiortho]

    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]

    if outdir is None:
        RESPATH = f"results/{rectPath}/{greedyPath}"
    else:
        RESPATH = outdir

    ### Initializa the nirb object
    nirb_off = nirbOffline(initCoarse=True, **config_nirb)
    nirb_off.initModel()

    start = time.time()

    Xi_train_path = RESPATH + "/sampling_train.sample"
    # if os.path.isfile(Xi_train_path):
    #     s = nirb_off.Dmu.sampling()
    #     N = s.readFromFile(Xi_train_path)
    #     Xi_train = s.getVector()
    #     if nirb_off.worldcomm.isMasterRank():
    #         print(f"[NIRB] Xi_train loaded from {Xi_train_path}")
    # else :
    #     Xi_train = generatedAndSaveSampling(nirb_off.Dmu, 200, path=Xi_train_path, samplingMode="log-random")

    Xi_train = generatedAndSaveSampling(nirb_off.Dmu, 200, path=Xi_train_path, samplingMode="log-random")

    nirb_off.generateOperators(coarse=True)

    if doGreedy:
        _, Xi_train, _ = nirb_off.initProblemGreedy(500, 1e-5, Nmax=config_nirb['nbSnapshots'], Xi_train=Xi_train, computeCoarse=True, samplingMode="random")
    else:
        Xi_train = nirb_off.initProblem(nbSnap, Xi_train=Xi_train)
    RIC = nirb_off.generateReducedBasis(tolerance=1.e-12)

    N = len(nirb_off.fineSnapShotList)

    tolortho = 1.e-8

    nirb_off.saveData(RESPATH, force=True)

    finish = time.time()

    print(f"proc {nirb_off.worldcomm.localRank()} Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized(tol=tolortho))
    print(f"proc {nirb_off.worldcomm.localRank()} Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    if convergence :
        Xi_test_path = RESPATH + "/sampling_test.sample"
        if os.path.isfile(Xi_test_path):
            s = nirb_off.Dmu.sampling()
            Ns = int(s.readFromFile(Xi_test_path))
            Xi_test = s.getVector()
            if nirb_off.worldcomm.isMasterRank():
                print(f"[NIRB] Xi_test loaded from {Xi_test_path}")
        else:
            Ns = 30
            Xi_test = generatedAndSaveSampling(nirb_off.Dmu, Ns, path=Xi_test_path, samplingMode="log-random")

        nirb_on = nirbOnline(**config_nirb)
        nirb_on.initModel()
        nirb_on.loadData(nbSnap=N, path=RESPATH)

        Err = nirb_off.checkConvergence(nirb_on, Ns=30, Xi_test=Xi_test)
        df = pd.DataFrame(Err)
        file = RESPATH + "/offlineError.csv"
        df.to_csv(file, index=False)

        if nirb_off.worldcomm.isMasterRank(): print(f"[NIRB] Offline error saved in {file}")


    perf = []
    perf.append(nbSnap)
    perf.append(finish-start)

    if doRectification:
        file = RESPATH + f'/nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}_rectif.dat'
    else :
        file = RESPATH + f'/nirbOffline_time_exec_np{nirb_off.worldcomm.globalSize()}.dat'
    WriteVecAppend(file, perf)

    info = nirb_off.getOfflineInfos()

    if nirb_off.worldcomm.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] doRectification : {nirb_off.doRectification}, doGreedy : {nirb_off.doGreedy}")
        print(f"[NIRB] Offline part Done !")

    fppc.Environment.saveTimers(display=True)