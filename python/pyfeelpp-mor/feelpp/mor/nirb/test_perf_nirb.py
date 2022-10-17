import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import time
import pandas as pd
from pathlib import Path
from nirb_perf import *

H = 0.5  # CoarseMeshSize
h = 0.1  # Fine mesh size
dim = 3
PWD = os.getcwd()
toolboxType='heat'
modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
modelsFolder = f"{PWD}/model/"
cfg_path = f"{modelsFolder}thermal-fin-3d/thermal-fin.cfg"
geo_path = f"{modelsFolder}thermal-fin-3d/fin.geo"
model_path = f"{modelsFolder}thermal-fin-3d/thermal-fin.json"
doRectification = False


def offline(N):

    start = time.time()
    nirb = nirbOffline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)
    nirb.initProblem(N)
    nirb.generateOperators()
    nirb.generateReducedBasis(regulParam=1.e-10)
    finish = time.time()

    nirb.saveData()

    perf = [N, finish-start]
    file = 'nirbOffline_time_exec.dat'
    WriteVecAppend(file, perf)

    print(f"[NIRB] Offline Elapsed time = ", finish-start)
    print(f"[NIRB] Offline part Done !")


def online_error_sampling():
    nirb = nirbOnline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)
    nirb.loadData()

    Nsample = 50
    errorN = ComputeErrorSampling(nirb, Nsample=Nsample)

    df = pd.DataFrame(errorN)
    if doRectification:
        file = Path(f"errorParamsRectif/errors{nirb.N}.csv")
    else:
        file = Path(f"errorParams/errors{nirb.N}.csv")
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


def online_time_measure():
    nirb = nirbOnline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)
    nirb.loadData()
    
    Dmu = nirb.Dmu
    Ns = 50
    s = Dmu.sampling()
    s.sampling(Ns, 'log-random')
    mus = s.getVector()

    time_toolbox_start = time.time()
    for mu in mus:
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
    time_toolbox_finish = time.time()
    time_toolbox = (time_toolbox_finish - time_toolbox_start) / Ns
    print(f"[NIRB online] Time to compute {Ns} solutions with toolbox = ", time_toolbox)

    time_nirb_start = time.time()
    for mu in mus:
        uHh = nirb.getOnlineSol(mu)
    time_nirb_finish = time.time()
    time_nirb = (time_nirb_finish - time_nirb_start) / Ns
    print(f"[NIRB online] Time to compute {Ns} solutions with NIRB = ", time_nirb)

    WriteVecAppend('nirbOnline_time_exec.dat', [nirb.N, time_toolbox, time_nirb])

    

if __name__ == '__main__':


    e = init_feelpp_environment(toolboxType, cfg_path)
    Ns = sys.argv[1:]

    for N in Ns:
        N = int(N)
        print("\n\n-----------------------------")
        print(f"[NIRB] Test with N = {N}")
        offline(N)
        online_error_sampling()
        online_time_measure()
