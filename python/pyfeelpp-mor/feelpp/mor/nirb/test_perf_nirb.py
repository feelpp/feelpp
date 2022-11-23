import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import pandas as pd
from pathlib import Path
from nirb_perf import *
import argparse

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def offline(config_nirb, RESPATH, doGreedy, N, Xi_train=None, regulParam=1e-10):
    """Run the offline phase of the NIRB method

    Args:
        config_nirb (dict): configuration of the NIRB
        RESPATH (str): path to the results
        doGreedy (bool): if True, the greedy algorithm is used to generate the reduced basis
        N (int): number of snapshots (limit number if greedy algo is used)
        Xi_train (list, optional): Set of parameters to train. Defaults to None.
        regulParam (_type_, optional): _description_. Defaults to 1e-10.
    """

    start = time.time()
    nirb = nirbOffline(**config_nirb, initCoarse=doGreedy)
    nirb.generateOperators(coarse=True)
    if doGreedy:
        nirb.initProblemGreedy(500, 1e-5, Nmax=N, computeCoarse=True, samplingMode="random")
    else:
        nirb.initProblem(N, Xi_train=Xi_train)
    nirb.generateReducedBasis(regulParam=regulParam)
    finish = time.time()

    nirb.saveData(RESPATH, force=True)

    perf = [N, finish-start]
    
    res = {}
    res['N'] = [N]
    res['nirb_offline'] = [finish-start]

    print(f"[NIRB offline] proc : {rank} Offline Elapsed time = ", res['nirb_offline'])

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] Offline part Done !")
    
    return res 


def online_error_sampling(nirb, RESPATH, Nsample=50, Xi_test=None, verbose=True, save=True):
    """Compute the error between the toolbox solution and the NIRB solution for a sampling of parameters
       Generates the file errorParams/errors${N}.csv containing the errors for each parameter, for a given value of N,
          corresponding to the size of the reduced basis.
    Args:
    -----
        nirb (nirbOnline): NIRB online object
        RESPATH (str): path to the results
        Nsample (int, optional): number of parameters to sample. Defaults to 50.
        Xi_test (list): list of parameters to test. Default is None, in which case the parameters are generated randomly
        verbose (bool, optional): if True, print the errors. Defaults to True.
        save (bool, optional): if True, save the errors in a csv file. Defaults to False.
    """

    Nsample = 50
    errorN = ComputeErrorSampling(nirb, Nsample=Nsample, Xi_test=Xi_test, h1=True)

    df = pd.DataFrame(errorN) 
    df['N'] = nirb.N 

    if save:
        if feelpp.Environment.isMasterRank():    
            file = Path(f"{RESPATH}/errors{Nsample}Params.csv")
            header = not os.path.isfile(file)
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


def online_time_measure(nirb, Nsample=50):
    """Measures the online time to compute solution, compared to the time taken by the toolbox
       Generates the file nirbOnline_time_exec.dat containing the time to compute the NIRB solution and the toolbox solution,
           according to the number of samples N
    
    Args:
    -----
        nirb_config (nirbOnline): nirb object
        Nsample (int, optional): number of parameters to sample. Defaults to 50.
    """
    Dmu = nirb.Dmu
    Ns = Nsample
    s = Dmu.sampling()
    s.sampling(Ns, 'log-random')
    mus = s.getVector()

    print("[NIRB online] Start toolbox time measure")
    time_toolbox_start = time.time()
    for mu in mus:
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
    time_toolbox_finish = time.time()
    time_toolbox = (time_toolbox_finish - time_toolbox_start) / Ns
    print(f"[NIRB online] Time to compute {Ns} solutions with toolbox = ", time_toolbox)

    print("[NIRB online] Start NIRB online time measure")
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

    model_path = config_nirb['model_path']
    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]
    RESPATH = f"results/{rectPath}/{greedyPath}"

    nbSnap=args.N
    if nbSnap==None:
        nbSnap = config_nirb['nbSnapshots']
        
    # Ns = sys.argv[1:]
    # Ns = args.N 
    Ns = [1, 2, 4, 6, 10, 12, 14, 16, 20, 25, 30, 35, 40, 45, 50] # number of basis functions 
    Ns = [1, 6]

    save = True 

    Dmu = loadParameterSpace(model_path)
    Xi_train = generatedAndSaveSampling(Dmu, 250, samplingMode="random")
    if False:
        Xi_test = generatedAndSaveSampling(Dmu, 50, path="Xi_test.sample", samplingMode="random")
    else:
        s = Dmu.sampling()
        N = s.readFromFile("Xi_test.sample")
        assert N == 50
        Xi_test = s.getVector()

    for N in Ns:
        N = int(N)
        print("\n\n-----------------------------")
        print(f"[NIRB] Test with N = {N}")
        resOffline = offline(config_nirb, RESPATH, doGreedy, N, Xi_train=None, regulParam=1e-10)
        nirb_on = nirbOnline(**config_nirb)
        err = nirb_on.loadData(path=RESPATH)
        assert err == 0, "loadData failed"
        online_error_sampling(nirb_on, RESPATH, Nsample=50, Xi_test=None, verbose=True, save=True)
        resOnline = online_time_measure(nirb_on, Nsample=50)

        if save:
            if feelpp.Environment.isMasterRank():
                res = resOnline
                res['nirb_offline'] = resOffline['nirb_offline']
                df = pd.DataFrame(res) 
                file = Path(f"{RESPATH}/nirb_time_exec.csv")
                header = not os.path.isfile(file)
                df.to_csv(file, mode='a', index=False, header=header)

    if feelpp.Environment.isMasterRank():
        print("=============================================")
        print(f"Run test_perf_nirb done, doRectification : {doRectification}, doGreedy : {doGreedy}")
        print(f"Result are stored in : {os.path.abspath(RESPATH)}")
        print("=============================================")
