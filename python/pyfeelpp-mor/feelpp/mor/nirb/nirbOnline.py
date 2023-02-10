from time import time
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import sys
import pandas as pd
from pathlib import Path
import os
from feelpp.mor.nirb.nirb_perf import *
import argparse


def run_online(config_nirb, path, nbSnap, Xi=None, Nb=None):
    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel()
    if Xi is None:
        s = nirb_on.Dmu.sampling()
        s.sampling(10, "random")
        Xi = s.getVector()

    if Nb is None: Nb = nbSnap

    err = nirb_on.loadData(path=path, nbSnap=nbSnap)
    assert err == 0, "Error while loading data"

    for mu in Xi:
        nirb_on.getOnlineSol(mu, Nb)

    return nirb_on


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=int, default=None)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)
    parser.add_argument("--greedy", help="With or without Greedy [default=0]", type=int, default=0)
    parser.add_argument("--exporter", help="Export nirb sol for vizualisation [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Get convergence error [default=1]", type=int, default=0)
    parser.add_argument("--compute-error", help="Compute errors of nirb solutions [default=0]", type=int, default=0)

    parser.add_argument("--Xi", help="Path to file containing the parameters", type=str, default=None)
    parser.add_argument("--Nparam", help="Number of parameters to test, if Xi is not given", type=int, default=10)

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']


    convergence = args.convergence != 0
    exporter = args.exporter != 0
    config_nirb['greedy-generation'] = args.greedy != 0
    computeError = args.compute_error != 0

    nbSnap = args.N
    if nbSnap == None:
        nbSnap = config_nirb['nbSnapshots']


    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]

    if outdir is None:
        RESPATH = f"results/{rectPath}/{greedyPath}"
    else:
        RESPATH = outdir

    start = time()
    nirb_on = run_online(config_nirb, path=RESPATH, nbSnap=nbSnap, Nb=nbSnap)
    finish = time()

    s = nirb_on.Dmu.sampling()
    if args.Xi is not None:
        N_train = s.readFromFile('./sampling.sample')
        Xi = s.getVector()
    else:
        s.sampling(args.Nparam, "random")
        Xi = s.getVector()



    if exporter:
        dirname = "nirbOnlineSolutions"
        nirb_on.initExporter(dirname, toolbox="fine")

    for i, mu in enumerate(Xi):
        print(f"[NIRB online] Getting online solution with mu = {mu}")
        uHh = nirb_on.getOnlineSol(mu)

        if computeError:
            uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
            error = nirb_on.normMat(uHh - uh)
            print(f"[NIRB] L2 norm between nirb online and toolbox sol = {error}")


        if exporter:
            nirb_on.exportField(uHh, f"uHhN_{i}")

    if exporter:
        nirb_on.saveExporter()

    perf = []
    perf.append(nirb_on.N)
    perf.append(finish-start)

    if doRectification:
        file = os.path.join(RESPATH, f'nirbOnline_time_exec_np{nirb_on.worldcomm.globalSize()}_rectif.dat')
    else:
        file = os.path.join(RESPATH, f'nirbOnline_time_exec_np{nirb_on.worldcomm.globalSize()}.dat')
    WriteVecAppend(file, perf)

    if convergence :
        Nsample = 50
        errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample, h1=True)
        # errorN = ComputeErrors(nirb_on, mu, h1=True)

        df = pd.DataFrame(errorN)
        df['N'] = nirb_on.N

        file = os.path.join(RESPATH, f"errors{Nsample}Params.csv")

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

    print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")
