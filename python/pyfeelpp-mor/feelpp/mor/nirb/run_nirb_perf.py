import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import init_feelpp_environment, generatedAndSaveSampling, merge_JsonFiles
import pandas as pd
from pathlib import Path
from nirb_perf import *
import argparse
import feelpp.core as fppc
from feelpp.core.timing import tic, toc

from mpi4py import MPI
comm = MPI.COMM_WORLD

def offline(nirb, RESPATH, doGreedy, N, Xi_train=None, regulParam=1e-10):
    """Run the offline phase of the NIRB method

    Args:
        nirb (class): offline nirb object
        RESPATH (str): path to the results
        doGreedy (bool): if True, the greedy algorithm is used to generate the reduced basis
        N (int): number of snapshots (limit number if greedy algo is used)
        Xi_train (list, optional): Set of parameters to train. Defaults to None.
        regulParam (_type_, optional): _description_. Defaults to 1e-10.
    """

    tic()
    nirb.generateOperators(coarse=True)
    if doGreedy:
        nirb.initProblemGreedy(500, 1e-3, Xi_train=Xi_train, Nmax=N, computeCoarse=True, samplingMode="log-random")
    else:
        nirb.initProblem(N, Xi_train=Xi_train)
    RIC = nirb.generateReducedBasis()

    nirb.saveData(RESPATH, force=True)

    RIC = np.array(RIC)
    file = "ric_offline.txt"
    np.savetxt(file, RIC)

    print(f"proc {nirb_off.worldcomm.localRank()} Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    offline_time = toc("NIRB: Offline time")

    comm.Barrier()

    res = {}
    res['N'] = [N]
    res['nirb_offline'] = [offline_time]

    if fppc.Environment.isMasterRank():
        print("[NIRB] Offline Elapsed time =", offline_time)
        print("[NIRB] Offline part Done !")

    return res


def online_error_sampling(nirb, RESPATH, Nb=None,  Nsample=50, Xi_test=None, verbose=True, save=True):
    """Compute the error between the toolbox solution and the NIRB solution for a sampling of parameters
       Generates the file errorParams/errors${N}.csv containing the errors for each parameter, for a given value of N,
          corresponding to the size of the reduced basis.
    Args:
    -----
        nirb (nirbOnline): NIRB online object
        RESPATH (str): path to the results
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used
        Nsample (int, optional): number of parameters to sample. Defaults to 50.
        Xi_test (list): list of parameters to test. Default is None, in which case the parameters are generated randomly
        verbose (bool, optional): if True, print the errors. Defaults to True.
        save (bool, optional): if True, save the errors in a csv file. Defaults to False.
    """

    errorN = ComputeErrorSampling(nirb, Nb=Nb, Nsample=Nsample, Xi_test=Xi_test, h1=True)

    df = pd.DataFrame(errorN)
    df['N'] = nirb.N

    if save:
        if fppc.Environment.isMasterRank():
            file = Path(f"{RESPATH}/errors{Nsample}Params.csv").absolute()
            header = not os.path.isfile(file)
            df.to_csv(file, mode='a', index=False, header=header)
            print(f"[NIRB] Convergence errors are saved in {file}")

    if verbose:
        if fppc.Environment.isMasterRank():
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


def online_time_measure(nirb, Nsample=50, Xi_test=None):
    """Measures the online time to compute solution, compared to the time taken by the toolbox
       Generates the file nirbOnline_time_exec.dat containing the time to compute the NIRB solution and the toolbox solution,
           according to the number of samples N

    Args:
    -----
        nirb_config (nirbOnline): nirb object
        Nsample (int, optional): number of parameters to sample. Defaults to 50.
    """
    if Xi_test is None :
        Dmu = nirb_on.Dmu
        s = Dmu.sampling()
        s.sample(Nsample, 'log-random')
        mus = s.getVector()
    else :
        mus = Xi_test

    toolbox_time = 0
    for mu in tqdm(mus, desc=f"[NIRB] Compute toolbox time mesure :", ascii=False, ncols=120):
        tic()
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
        toolbox_time += toc("NIRB: Toolbox time")
    time_toolbox = toolbox_time / Nsample

    nirb_time = 0
    for mu in tqdm(mus,desc=f"[NIRB] Compute online time mesure :", ascii=False, ncols=120):
        tic()
        uHh = nirb.getOnlineSol(mu)
        nirb_time += toc("NIRB: Online time")
    time_nirb = nirb_time / Nsample

    comm.Barrier()

    res = {}
    res['N'] = [nirb.N]
    res['nirb_online'] = [time_nirb]
    res['toolbox'] = [time_toolbox]

    if fppc.Environment.isMasterRank():
        print(f"[NIRB online] Time to compute {Nsample} solutions with NIRB = ", time_nirb)
        print(f"[NIRB online] Time to compute {Nsample} solutions with toolbox = ", time_toolbox)

    return res



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--Ntest", help="Number of sampling test parameter [default=10]", type=int, default=10)
    parser.add_argument("--Nbase", help="Max number of basis function [default=100]", type=int, default=100)
    parser.add_argument("--greedy", help="With or without Greedy [default=0]", type=int, default=0)
    parser.add_argument("--savetime", help="Save or not the execution time [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Get convergence error [default=1]", type=int, default=1)
    parser.add_argument("--regul", help="Regularization parameter [default=1e-10]", type=float, default=1e-10)

    ## get parser args
    args = parser.parse_args()
    config_file = args.config_file

    greedy = args.greedy
    save_time = args.savetime
    conv = args.convergence
    Nsample = args.Ntest
    Nbase = args.Nbase
    lambda_regul = args.regul

    bo = [False, True]
    timeExec = bo[save_time]
    convergence = bo[conv]

    ## Init feelpp
    cfg = fppc.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType'][0]
    e = init_feelpp_environment(toolboxType, config_file)
    nirb_file = fppc.Environment.expand(cfg['nirb']['filename'][0])
    config_nirb = fppc.readJson(nirb_file)['nirb']

    ## Get model and data path
    config_nirb['greedy-generation'] = bo[greedy]
    model_path = config_nirb['model_path']
    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]
    RESPATH = f"results/{rectPath}/{greedyPath}"

    size = fppc.Environment.numberOfProcessors()

    ### For time mesure in // computing

    ## square4 2D :
    # config_nirb['coarsemesh_path'] = f"$cfgdir/meshFiles/squareCoarse_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/meshFiles/squareFine_p{size}.json"
    Xi_train_path = RESPATH + f"/sampling_train4_N200.sample"
    Xi_test_path = RESPATH + f"/sampling_test4_Ns{Nsample}.sample"

    ## square9 2D
    # config_nirb['coarsemesh_path'] = f"$cfgdir/square9mesh/squareCoarse_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/square9mesh/squareFine_p{size}.json"
    # Xi_train_path = RESPATH + f"/sampling_train9_N200.sample"
    # Xi_test_path = RESPATH + f"/sampling_test9_Ns{Nsample}.sample"

    ## thermal fin 3D
    # config_nirb['coarsemesh_path'] = f"$cfgdir/meshFiles/coarseMesh_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/meshFiles/fineMesh_p{size}.json"
    # Xi_train_path = RESPATH + f"/sampling_train3dfin_N200.sample"
    # Xi_test_path = RESPATH + f"/sampling_test3dfin_Ns{Nsample}.sample"

    pas = 4
    baseList = range(1, Nbase, pas)

    Dmu = loadParameterSpace(model_path)
    # Xi_train_path = RESPATH + f"/sampling_train4_N200.sample"
    # Xi_test_path = RESPATH + f"/sampling_test4_Ns{Nsample}.sample"

    if os.path.isfile(Xi_train_path):
        s = Dmu.sampling()
        N = s.readFromFile(Xi_train_path)
        Xi_train = s.getVector()
        if fppc.Environment.isMasterRank():
            print(f"[NIRB] Xi_train loaded from {Xi_train_path}")
    else:
        Xi_train = generatedAndSaveSampling(Dmu, 200, path=Xi_train_path, samplingMode="log-random")

    if os.path.isfile(Xi_test_path):
        s = Dmu.sampling()
        N = s.readFromFile(Xi_test_path)
        assert N == Nsample, f"Given size of sampling test {Nsample} # loaded sampling size {N}"
        Xi_test = s.getVector()
        if fppc.Environment.isMasterRank():
            print(f"[NIRB] Xi_test loaded from {Xi_test_path}")
    else :
        Xi_test = generatedAndSaveSampling(Dmu, Nsample, path=Xi_test_path, samplingMode="log-random")

    ## generate nirb offline and online object :
    nirb_off = nirbOffline(**config_nirb, initCoarse=doGreedy)
    full_model = merge_JsonFiles(cfg[toolboxType]["json.filename"])
    nirb_off.initModel(model=full_model)
    resOffline = offline(nirb_off, RESPATH, doGreedy, baseList[-1], Xi_train=Xi_train, regulParam=lambda_regul)
    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel(model=full_model)

    # Nglob = nirb_off.N
    # err = nirb_on.loadData(nbSnap=Nglob, path=RESPATH)
    # assert err == 0, "loadData failed"



    comm.Barrier()

    if convergence :
        ## Get convergence error
        for N in baseList:
            N = int(N)
            if fppc.Environment.isMasterRank():
                print("\n\n-----------------------------")
                print(f"[NIRB] Test with N = {N}")
            err = nirb_on.loadData(nbSnap=N, path=RESPATH)
            assert err == 0, "loadData failed"
            online_error_sampling(nirb_on, RESPATH, Nsample=Nsample, Xi_test=Xi_test, verbose=False, save=True)

    if timeExec:
        ## Get online time mesure
        err = nirb_on.loadData(nbSnap=baseList[-1], path=RESPATH)
        assert err == 0, "loadData failed"
        resOnline = online_time_measure(nirb_on, Nsample=Nsample)
        if fppc.Environment.isMasterRank():
            res = resOnline
            res['nirb_offline'] = resOffline['nirb_offline']
            res['nproc'] = size
            df = pd.DataFrame(res)
            # save file
            parent3 = list(Path(RESPATH).absolute().parents)[3]
            file = Path(f"{parent3}/nirb_parallel_time_exec.csv")
            header = not os.path.isfile(file)
            df.to_csv(file, mode='a', index=False, header=header)
            print(f"[NIRB] Execution time mesure are saved in {file}")

    if fppc.Environment.isMasterRank():
        print("=============================================")
        print(f"Run test_perf_nirb done, doRectification : {doRectification}, doGreedy : {doGreedy}")
        print(f"Result are stored in : {os.path.abspath(RESPATH)}")
        print(f"Set convergence study = {convergence}")
        print(f"Set execution time = {timeExec}")
        print(f"Elapsed time offline = {resOffline}")
        print("=============================================")

    fppc.Environment.saveTimers(display=True)