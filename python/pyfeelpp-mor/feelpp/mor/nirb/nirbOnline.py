from time import time
import feelpp.core as fppc
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import sys
import pandas as pd
from pathlib import Path
import os
from feelpp.mor.nirb.nirb_perf import *
import argparse
from feelpp.core.timing import tic, toc

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=None)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)
    parser.add_argument("--greedy", help="With or without Greedy [default=0]", type=int, default=0)
    parser.add_argument("--exporter", help="Export nirb sol for vizualisation [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Get convergence error [default=1]", type=int, default=0)

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = fppc.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = fppc.Environment.expand(cfg['nirb']['filename'])
    config_nirb = fppc.readJson(nirb_file)['nirb']


    greedy = args.greedy
    expo = args.exporter
    conv = args.convergence

    bo = [False, True]
    exporter = bo[expo]
    convergence = bo[conv]

    nbSnap=args.N
    if nbSnap==None:
        nbSnap = config_nirb['nbSnapshots']

    config_nirb['greedy-generation'] = bo[greedy]

    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]

    if outdir is None:
        RESPATH = f"results/{rectPath}/{greedyPath}"
    else:
        RESPATH = outdir

    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel()

    start= time()

    mu = nirb_on.Dmu.element()
    err = nirb_on.loadData(path=RESPATH, nbSnap=nbSnap)
    assert err == 0, "Error while loading data"
    tic()
    uHh_r = nirb_on.getOnlineSol(mu, doRectification=True)
    toc("NIRB w/ rectification")
    tic()
    uHh = nirb_on.getOnlineSol(mu, doRectification=False)
    toc("NIRB w/o rectification")
    finish = time()

    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    error_r = nirb_on.normMat(uHh_r - uh)
    error = nirb_on.normMat(uHh - uh)

    if nirb_on.worldcomm.isMasterRank():
        print(f"[NIRB] L2 norm between nirb online rectified and toolbox sol = {error_r}")
        print(f"[NIRB] L2 norm between nirb online and toolbox sol = {error}")

    if exporter:
        if nirb_on.worldcomm.isMasterRank():
            print(f"[NIRB] Exporting nirb sol for vizualisation")
        dirname = "nirbSol"
        nirb_on.initExporter(dirname, toolbox="fine")
        nirb_on.exportField(uh, "uh")
        nirb_on.exportField(uHh_r, "uNirb_r")
        nirb_on.exportField(uHh, "uNirb")
        nirb_on.saveExporter()



    perf = []
    perf.append(nirb_on.N)
    perf.append(finish-start)

    if doRectification:
        file=RESPATH+f'/nirbOnline_time_exec_np{nirb_on.worldcomm.globalSize()}_rectif.dat'
    else :
        file=RESPATH+f'/nirbOnline_time_exec_np{nirb_on.worldcomm.globalSize()}.dat'
    WriteVecAppend(file,perf)


    if convergence :
        Nsample = 50
        errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample, h1=True)
        # errorN = ComputeErrors(nirb_on, mu, h1=True)

        df = pd.DataFrame(errorN)
        df['N'] = nirb_on.N

        file =RESPATH +f"/errors{Nsample}Params.csv"

        header = not os.path.isfile(file)
        df.to_csv(file, mode='a', index=False, header=header)

        if nirb_on.worldcomm.isMasterRank():
            print(f"[NIRB online] with {nirb_on.N} snapshots ")
            print(f"[NIRB online] computed errors for {df.shape[0]} parameters ")
            data_mean = df.mean(axis=0)
            print("[NIRB online] Mean of errors ")
            print(data_mean)
            data_min = df.min(axis=0)
            print("[NIRB online] Min of errors ")
            print(data_min)
            data_max = df.max(axis=0)
            print("[NIRB online] Max of errors ")
            print(data_max)


    if nirb_on.worldcomm.isMasterRank():
        print(f"[NIRB] Online Elapsed time =", finish-start)
        print(f"[NIRB] Online part Done !!")

    fppc.Environment.saveTimers(display=True)