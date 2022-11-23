from time import time
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import sys
import pandas as pd
from pathlib import Path
import os 
from feelpp.mor.nirb.nirb_perf import *
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='NIRB Online')
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
    toolboxType = config_nirb['toolboxType']

    nbSnap=args.N
    if nbSnap==None:
        nbSnap = config_nirb['nbSnapshots']

    if outdir is None:
        RESPATH = config_nirb["outdir"]
    else:
        RESPATH = outdir

    start= time()
    
    nirb_on = nirbOnline(**config_nirb)

    mu = nirb_on.Dmu.mumin()
    err = nirb_on.loadData(path=RESPATH, nbSnap=nbSnap)
    assert err == 0, "Error while loading data"
    uHh = nirb_on.getOnlineSol(mu)

    finish = time()
    
    perf = []
    perf.append(nirb_on.N)
    perf.append(finish-start)

    doRectification = config_nirb['doRectification']
    Nsample = 50

    if doRectification:
        file='nirbOnline_time_exec_rectif.txt'
    else :
        file='nirbOnline_time_exec.txt'
    WriteVecAppend(file,perf)

    # error1 = ComputeErrors(nirb_on, mu)
    errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample, h1=True)
    # errorR = ComputeErrors(nirb_on, mu)

    df = pd.DataFrame(errorN)
    df['N'] = nirb_on.N

    file =f"errors{Nsample}Params.csv"

    # file = f'errorRelative.csv'

    header = True
    if os.path.isfile(file):
        header=False
    df.to_csv(file, mode='a', index=False, header=header)

    # if feelpp.Environment.isMasterRank():
    #     print(f"[NIRB online] with {nirb_on.N} snapshots ")
    #     print(f"[NIRB online] computed errors for {df.shape[0]} parameters ")
    #     data_mean = df.mean(axis=0)
    #     print("[NIRB online] Mean of errors ")
    #     print(data_mean)
    #     data_min = df.min(axis=0)
    #     print("[NIRB online] Min of errors ")
    #     print(data_min)
    #     data_max = df.max(axis=0)
    #     print("[NIRB online] Max of errors ")
    #     print(data_max)
       


    print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")