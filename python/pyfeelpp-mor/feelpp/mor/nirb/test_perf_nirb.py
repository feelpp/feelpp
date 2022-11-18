import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import time
import pandas as pd
from pathlib import Path
from feelpp.mor.nirb.nirb_perf import *
import argparse

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def offline(config_nirb, N, regulParam=1.e-10):

    start = time.time()
    nirb = nirbOffline(**config_nirb)
    nirb.initProblem(N)
    nirb.generateOperators()
    nirb.generateReducedBasis(regulParam=regulParam)
    finish = time.time()
    nirb.saveData(force=True)

    perf = [N, finish-start]

    res = {}
    res['N'] = [N]
    res['nirb_offline'] = [finish-start]

    print(f"[NIRB offline] proc : {rank} Offline Elapsed time = ", res['nirb_offline'])

    if feelpp.Environment.isMasterRank():
        file = 'nirbOffline_time_exec.dat'
        WriteVecAppend(file, perf)
        # print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] Offline part Done !")
    
    return res 

def online_error_sampling(config_file, nbSnap, Nsample=50, verbose=True, save=False):

    nirb = nirbOnline(**config_file)
    nirb.loadData(nbSnap=nbSnap)

    errorN = ComputeErrorSampling(nirb, Nsample=Nsample, h1=True)

    df = pd.DataFrame(errorN)
    df['N'] = nirb.N 

    if save:
        if feelpp.Environment.isMasterRank():    
            file = Path(f"errors{Nsample}Params.csv")
            header = True
            if os.path.isfile(file):
                header=False
            df.to_csv(file, mode='a', index=False, header=header)
            print(f"[NIRB online] Dataframe saved in {file}")
        
    if verbose:
        if feelpp.Environment.isMasterRank():
            print("[NIRB online] all computed errors ")
            data_mean = df.mean(axis=0)
            print("[NIRB online] Mean of errors ")
            print(data_mean)
            data_min = df.min(axis=0)
            print("[NIRB online] Min of errors ")
            print(data_min)
            data_max = df.max(axis=0)
            print("[NIRB online] Max of errors ")
            print(data_max)


def online_time_measure(config_nirb, nbSnap, Nsample=50):
    
    nirb = nirbOnline(**config_nirb)
    nirb.loadData(nbSnap=nbSnap)

    Dmu = nirb.Dmu
    Ns = Nsample
    s = Dmu.sampling()
    s.sampling(Ns, 'log-random')
    mus = s.getVector()

    time_toolbox_start = time.time()
    for mu in mus:
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
    time_toolbox_finish = time.time()
    time_toolbox = (time_toolbox_finish - time_toolbox_start) / Ns

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB online] Time to compute {Ns} solutions with toolbox = ", time_toolbox)

    time_nirb_start = time.time()
    for mu in mus:
        uHh = nirb.getOnlineSol(mu)
    time_nirb_finish = time.time()
    time_nirb = (time_nirb_finish - time_nirb_start) / Ns

    res = {}
    res['N'] = [nirb.N]
    res['nirb_online'] = [time_nirb]
    res['toolbox'] = [time_toolbox]

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB online] Time to compute {Ns} solutions with NIRB = ", time_nirb)
        WriteVecAppend('nirbOnline_time_exec.dat', [nirb.N, time_toolbox, time_nirb])

    return res 
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--N", help="Number of initial snapshots [default=10]", type=list, default=None)

    args = parser.parse_args()
    config_file = args.config_file

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']
    toolboxType = config_nirb['toolboxType']

    nbSnap=args.N
    if nbSnap==None:
        nbSnap = config_nirb['nbSnapshots']
        
    # Ns = sys.argv[1:]
    # Ns = args.N 
    Ns = [1, 2, 4, 6, 10, 12, 14, 16, 20, 25, 30, 35, 40, 45, 50] # number of basis functions 
    Ns = [1, 2]

    save = True 

    for N in Ns:
        N = int(N)
        print("\n\n-----------------------------")
        print(f"[NIRB] Test with N = {N}")
        resOffline = offline(config_nirb, N, 1.e-10)
        online_error_sampling(config_nirb, N, Nsample=30, verbose=True)
        resOnline = online_time_measure(config_nirb, N, Nsample=30)

        if save:
            if feelpp.Environment.isMasterRank():
                res = resOnline
                res['nirb_offline'] = resOffline['nirb_offline']
                df = pd.DataFrame(res) 
                file = Path(f"nirb_time_exec.csv")
                header = True
                if os.path.isfile(file):
                    header=False
                df.to_csv(file, mode='a', index=False, header=header)

