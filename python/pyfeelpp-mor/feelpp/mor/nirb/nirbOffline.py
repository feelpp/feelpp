import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import json
import argparse
from test_perf_nirb import online_error_sampling


if __name__ == "__main__":
    PWD = os.getcwd()

    parser = argparse.ArgumentParser(description='NIRB Offline')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=None)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']

    nbSnap=args.N
    if nbSnap==None:
        nbSnap = config_nirb['nbSnapshots']

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

    start = time.time()
    ###
    # Only once: generate and save a sampling
    # if False:
    #     generatedAndSaveSampling(nirb_off.Dmu, 100)
    #     sys.exit(0)
    ###
    # If wanted: load the savec sampling to use it in algorithm generation
    Xi_train = None
    # if True:
    #     s = nirb_off.Dmu.sampling()
    #     N_train = s.readFromFile('./sampling.sample')
    #     Xi_train = s.getVector()

    nirb_off.generateOperators(coarse=True)

    if doGreedy:
        _, Xi_train, _ = nirb_off.initProblemGreedy(100, 1e-5, Nmax=config_nirb['nbSnapshots'], Xi_train=Xi_train, computeCoarse=True, samplingMode="random")
    else:
        Xi_train = nirb_off.initProblem(nbSnap)
    nirb_off.generateReducedBasis(regulParam=1.e-10)
    # nirb_off.orthonormalizeL2()
    # nirb_off.orthonormalizeH1()
    nirb_off.saveData(RESPATH, force=True)

    finish = time.time()

    # if feelpp.Environment.isMasterRank():
    #     print("Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    #     print("Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    perf = []
    perf.append(nbSnap)
    perf.append(finish-start)

    if doRectification:
        file=RESPATH+f'/nirbOffline_time_exec_np{size}_rectif.dat'
    else :
        file=RESPATH+f'/nirbOffline_time_exec_np{size}.dat'
    WriteVecAppend(file, perf)

    info = nirb_off.getInformations()

    if feelpp.Environment.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] Offline part Done !")
