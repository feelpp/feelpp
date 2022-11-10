import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import pandas as pd
from pathlib import Path
from nirb_perf import *



def offline(nirb_config, RESPATH, doGreedy, N, Xi_train=None):
    """Generated the reduced basis for the given N
       Generated a file nirbOffline_time_exec.dat containing the time to compute the reduced basis,
         according to the number of samples N

    Args:
    -----
        nirb_config (dict): configuration of the NIRB
        RESPATH (str): path to the results
        doGreedy (bool): if True, the greedy algorithm is used to generate the reduced basis
        N (int): number of snapshots (limit number if greedy algo is used)

    Return:
    --------
        nirb (nirbOnline): NIRB object
    """

    start = time.time()
    nirb = nirbOffline(**nirb_config, initCoarse=doGreedy)

    nirb.generateOperators(coarse=True)
    if doGreedy:
        nirb.initProblemGreedy(500, 1e-5, Nmax=N, computeCoarse=True, samplingMode="random")
    else:
        nirb.initProblem(N, Xi_train=Xi_train)
    nirb.generateReducedBasis(regulParam=1.e-10)
    finish = time.time()

    nirb.saveData(RESPATH, force=True)

    perf = [nirb.N, finish-start]
    file = RESPATH + '/nirbOffline_time_exec.dat'
    WriteVecAppend(file, perf)

    print(f"[NIRB] Offline Elapsed time = ", finish-start)
    print(f"[NIRB] Offline part Done !")
    return nirb


def online_error_sampling(nirb, RESPATH, Xi_test=None):
    """Compute the error between the toolbox solution and the NIRB solution for a sampling of parameters
       Generates the file errorParams/errors${N}.csv containing the errors for each parameter, for a given value of N,
          corresponding to the size of the reduced basis.

    Args:
    -----
        nirb_config (nirbOnline): configuration of the NIRB
        RESPATH (str): path to the results
        Xi_test (list): list of parameters to test. Default is None, in which case the parameters are generated randomly
    """

    Nsample = 50
    errorN = ComputeErrorSampling(nirb, Nsample=Nsample, Xi_test=Xi_test)

    df = pd.DataFrame(errorN)
    file = Path(f"{RESPATH}/errorParams/errors{nirb.N}.csv")
    file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(file, index=False)
    print(f"[NIRB online] Dataframe saved in {file}")

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


def online_time_measure(nirb, RESPATH):
    """Measures the online time to compute solution, compared to the time taken by the toolbox
       Generates the file nirbOnline_time_exec.dat containing the time to compute the NIRB solution and the toolbox solution,
           according to the number of samples N
    
    Args:
    -----
        nirb_config (nirbOnline): nirb object
        RESPATH (str): path to the results
    """
    Dmu = nirb.Dmu
    Ns = 50
    s = Dmu.sampling()
    s.sampling(Ns, 'log-random')
    mus = s.getVector()

    print("[NIRB online] Start toolbox time measure")
    time_toolbox_start = time.time()
    for mu in mus:
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
    time_toolbox_finish = time.time()
    time_toolbox = (time_toolbox_finish - time_toolbox_start) / Ns
    print(f"[NIRB online] Average time to compute {Ns} solutions with toolbox = ", time_toolbox)

    print("[NIRB online] Start NIRB online time measure")
    time_nirb_start = time.time()
    for mu in mus:
        uHh = nirb.getOnlineSol(mu)
    time_nirb_finish = time.time()
    time_nirb = (time_nirb_finish - time_nirb_start) / Ns
    print(f"[NIRB online] Average time to compute {Ns} solutions with NIRB = ", time_nirb)

    WriteVecAppend(RESPATH+'/nirbOnline_time_exec.dat', [nirb.N, time_toolbox, time_nirb])


if __name__ == '__main__':

    PWD = os.getcwd()
    toolboxType = "heat"
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}thermal-fin-3d/thermal-fin.cfg"
    geo_path = f"{modelsFolder}thermal-fin-3d/fin.geo"
    model_path = f"{modelsFolder}thermal-fin-3d/thermal-fin.json"
    doRectification = True
    doGreedy = False
    nirb_config = feelpp.readJson(model_path)['nirb']
    nirb_config['doRectification'] = doRectification
    nirb_config['doGreedy'] = doGreedy

    r = ["noRect","Rect"][doRectification]
    g = ["noGreedy","Greedy"][doGreedy]
    RESPATH = f"RESULTS/{r}/{g}"

    e = init_feelpp_environment(toolboxType, cfg_path)
    Ns = sys.argv[1:]

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
        nirb = offline(nirb_config, RESPATH, doGreedy, N, Xi_train=Xi_train)
        nirb_on = nirbOnline(**nirb_config)
        err = nirb_on.loadData(path=RESPATH)
        assert err == 0, "loadData failed"
        online_error_sampling(nirb_on, RESPATH, Xi_test=Xi_test)
        online_time_measure(nirb_on, RESPATH)

    if feelpp.Environment.isMasterRank():
        print("=============================================")
        print(f"Run test_perf_nirb done, doRectification : {doRectification}, doGreedy : {doGreedy}")
        print(f"Result are stored in : {RESPATH}")
        print("=============================================")