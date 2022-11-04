import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import json
import argparse
from test_perf_nirb import online_error_sampling, RESPATH


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='NIRB Offline')
    parser.add_argument('--config-file', type=str, help='path to cfg file')

    args = parser.parse_args()
    config_file = args.config_file

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']

    doRectification  = config_nirb['doRectification'] 
    nbSnap = config_nirb['nbSnapshots']

    start = time.time()

    nirb_off = nirbOffline(initCoarse=True, **config_nirb)
    
    ### 
    # Only once: generate and save a sampling
    if False:
        generatedAndSaveSampling(nirb_off.Dmu, 100)
        sys.exit(0)
    ###
    # If wanted: load the savec sampling to use it in algorithm generation
    Xi_train = None
    if True:
        s = nirb_off.Dmu.sampling()
        N_train = s.readFromFile('./sampling.sample')
        Xi_train = s.getVector()

    nirb_off.generateOperators(coarse=True)

    doGreedy = config_nirb['greedy-generation']
    if doGreedy:
        res = nirb_off.initProblemGreedy(100, 1e-5, Nmax=6, Xi_train=Xi_train, computeCoarse=True, samplingMode="random")
    else:
        nirb_off.initProblem(nbSnap)
    nirb_off.generateReducedBasis(regulParam=1.e-10)
    nirb_off.saveData(force=True)

    
    finish = time.time()

    # nirb_off.orthonormalizeL2()
    # nirb_off.orthonormalizeH1()
    nirb_off.saveData(RESPATH)

    online_error_sampling(res[1])

    print("Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    print("Is H1 orthonormalized ? ", nirb_off.checkH1Orthonormalized())

    perf = []
    perf.append(nbSnap)
    perf.append(finish-start)

    file='nirbOffline_time_exec.txt'
    WriteVecAppend(file,perf)

    if doRectification:
        file='nirbOffline_time_exec_rectif.txt'
    else :
        file='nirbOffline_time_exec.txt'
    WriteVecAppend(file,perf)

    info = nirb_off.getInformations()

    if feelpp.Environment.isMasterRank():
        print(json.dumps(info, sort_keys=True, indent=4))
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] Offline part Done !")
