import sys
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment, generatedAndSaveSampling
import time
import pandas as pd
from pathlib import Path
from nirb_perf import *
import argparse
from os.path import dirname, basename, isfile, join

from mpi4py import MPI
comm = MPI.COMM_WORLD

def offline(nirb, RESPATH, doGreedy, N, Xi_train=None):
    """Run the offline phase of the NIRB method

    Args:
        nirb (class): offline nirb object
        RESPATH (str): path to the results
        doGreedy (bool): if True, the greedy algorithm is used to generate the reduced basis
        N (int): number of snapshots (limit number if greedy algo is used)
        Xi_train (list, optional): Set of parameters to train. Defaults to None.
    """

    start = time.time()
    nirb.generateOperators(coarse=True)
    if doGreedy:
        nirb.initProblemGreedy(500, 1e-3, Xi_train=Xi_train, Nmax=N, computeCoarse=True, samplingMode="log-random")
    else:
        nirb.initProblem(N, Xi_train=Xi_train)
    RIC = nirb.generateReducedBasis()

    nirb.saveData(RESPATH, force=True)
    
    # RIC = np.array(RIC)
    # file = "ric_offline.txt"
    # np.savetxt(file,RIC)

    print(f"proc {nirb_off.worldcomm.localRank()} Is L2 orthonormalized ?", nirb_off.checkL2Orthonormalized())
    finish = time.time()
    
    comm.Barrier()

    res = {}
    res['N'] = [N]
    res['nirb_offline'] = [finish-start]

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB] Offline Elapsed time = ", finish-start)
        print(f"[NIRB] Offline part Done !")
    
    return res 


def online_error_sampling(nirb, RESPATH, errorfile=None, Nb=None,  Nsample=50, Xi_test=None, verbose=True, save=True, regulParam=1.e-10):
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
        regulParam(float, optionnal) : regularization parameter for the rectification preprocess. Defaults to 1.e-10
    """
    
    errorN = ComputeErrorSampling(nirb, Nb=Nb, Nsample=Nsample, Xi_test=Xi_test, h1=True, regulParam=regulParam)

    df = pd.DataFrame(errorN) 
    df['N'] = Nb

    if save:
        if feelpp.Environment.isMasterRank():
            if not os.path.isdir(RESPATH):
                os.makedirs(RESPATH)

            if errorfile is None :   
                file = Path(f"{RESPATH}/errors{Nsample}Params.csv").absolute()
            else : 
                if not RESPATH[-1] == "/":
                    file = RESPATH +'/'+ errorfile 
                else :
                    file = RESPATH + errorfile 

            header = not os.path.isfile(file)
            df.to_csv(file, mode='a', index=False, header=header)
            print(f"[NIRB] Convergence errors are saved in {file}")
        
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
        s.sampling(Nsample, 'log-random')
        mus = s.getVector()
    else :
        mus = Xi_test


    time_toolbox_start = time.time()
    for mu in tqdm(mus,desc=f"[NIRB] Compute toolbox time mesure :", ascii=False, ncols=120):
        uh = nirb.getToolboxSolution(nirb.tbFine, mu)
    time_toolbox_finish = time.time()
    time_toolbox = (time_toolbox_finish - time_toolbox_start) / Nsample

    time_nirb_start = time.time()
    for mu in tqdm(mus,desc=f"[NIRB] Compute online time mesure :", ascii=False, ncols=120):
        uHh = nirb.getOnlineSol(mu)
    time_nirb_finish = time.time()
    time_nirb = (time_nirb_finish - time_nirb_start) / Nsample

    comm.Barrier()

    res = {}
    res['N'] = [nirb.N]
    res['nirb_online'] = [time_nirb]
    res['toolbox'] = [time_toolbox]

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB online] Time to compute {Nsample} solutions with NIRB = ", time_nirb)
        print(f"[NIRB online] Time to compute {Nsample} solutions with toolbox = ", time_toolbox)

    return res

def loadSampling(Dmu, path='./', idmodel='s4', Ntest=50, Ntrain=200, samplingMode='log-random'):
    """upload sampling parameter for trainning and testing 

    Parameters
    ----------
    path : str
        path to load or save sampling 
    idmodel : str, optional
        identifiant of the model, by default 's4'
    Ntest : int, optional
        number of test parameter, by default 50
    Ntrain : int, optional
        number of train parameter, by default 200
    """

    Xtrain_path = path + f"/sampling_train_{idmodel}_N{Ntrain}.sample"
    Xtest_path = path + f"/sampling_test_{idmodel}_N{Ntest}.sample"

    
    if os.path.isfile(Xtest_path) and os.path.isfile(Xtrain_path) :
        s = Dmu.sampling()
        N = s.readFromFile(Xtest_path)
        assert N==Ntest, f"Given size of sampling test {Ntest} # loaded sampling size {N}"
        Xi_test = s.getVector() 
        N = s.readFromFile(Xtrain_path)
        assert N==Ntrain, f"Given size of sampling train {Ntrain} # loaded sampling size {N}"
        Xi_train = s.getVector() 
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB] Sampling loaded from path : {path}")  
    else :
        Xi_train, Xi_test = SamplingPreProcess(Dmu,Ntrain=Ntrain, Ntest=Ntest, path=path, idmodel=idmodel, samplingMode=samplingMode)

    return Xi_train, Xi_test
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--Ntest", help="Number of sampling test parameter [default=10]", type=int, default=10)
    parser.add_argument("--Nbase", help="Max number of basis function [default=100]", type=int, default=100)
    parser.add_argument("--greedy", help="Wether to get Greedy or not [default=0]", type=int, default=0)
    parser.add_argument("--timeexec", help="Wether to get the execution time or not [default=0]", type=int, default=0)
    parser.add_argument("--convergence", help="Wther to get convergence error [default=1]", type=int, default=1)
    parser.add_argument("--idmodel", help="identifiant of the model [default='s4'(for square 4)]", type=str, default="s4")
    parser.add_argument("--regulparam", help="Wether to get conv error in respect of regularization parameter [default=0]", type=int, default=0)

    ## get parser args 
    args = parser.parse_args()
    config_file = args.config_file 
    Nsample = args.Ntest
    Nbase = args.Nbase

    bo = [False, True]
    timeExec = bo[args.timeexec]
    convergence = bo[args.convergence]
    regulParameter = bo[args.regulparam]

    ## Init feelpp 
    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)
    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']

    ## Get model and data path 
    config_nirb['greedy-generation'] = bo[args.greedy]
    model_path = config_nirb['model_path']
    doGreedy = config_nirb['greedy-generation']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]
    RESPATH = f"results/{rectPath}/{greedyPath}"
    if feelpp.Environment.isMasterRank():
        if not os.path.exists(RESPATH):
            os.makedirs(RESPATH)
    
    size = feelpp.Environment.numberOfProcessors()

    
    ## square4 2D :
    # config_nirb['coarsemesh_path'] = f"$cfgdir/meshFiles/squareCoarse_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/meshFiles/squareFine_p{size}.json"
    # idmodel = 's4'
    
    ## square9 2D 
    # config_nirb['coarsemesh_path'] = f"$cfgdir/square9mesh/squareCoarse_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/square9mesh/squareFine_p{size}.json"
    # idmodel ='s9'
    
    ## thermal fin 3D 
    # config_nirb['coarsemesh_path'] = f"$cfgdir/meshFiles/coarseMesh_p{size}.json"
    # config_nirb['finemesh_path'] = f"$cfgdir/meshFiles/fineMesh_p{size}.json"
    

    pas = 3
    baseList = range(10,Nbase, pas)

    Dmu = loadParameterSpace(model_path)
    ## get distinct training and testing sample 
    idmodel = args.idmodel 
    Xi_train, Xi_test = loadSampling(Dmu, path=RESPATH, idmodel=idmodel, Ntest=Nsample)

    ## generate nirb offline and online object :  
    # nirb_off = nirbOffline(**config_nirb, initCoarse=doGreedy)
    # nirb_off.initModel()
    # resOffline = offline(nirb_off, RESPATH, doGreedy, baseList[-1], Xi_train=Xi_train)
    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel()

    ## load offline datas only once 
    lmd = 1.e-10
    # Nglob = nirb_off.N
    Nglob = 97
    err = nirb_on.loadData(nbSnap=Nglob, path=RESPATH, regulParam=lmd)
    assert err == 0, "loadData failed"
    idd = str(idmodel) + f"lmd{10}"
    errorfile = f"errors{Nsample}Params_{idd}.csv"

    comm.Barrier()
    Nl = [10]
    ## for n in Nl, Lambda = 1.E-n
    if convergence : 
        if regulParameter : 
            ## test convergence i respect of regularization parameter 
            Nl = [0, 1, 3, 5, 7, 9, 10, 11, 12, 17]
            del nirb_on.RectificationMat[Nglob] 
            for lm in Nl :
                lmd = 1./10**lm 
                idd = str(idmodel) + f"lmd{lm}"
                errorfile = f"errors{Nsample}Params_{idd}.csv"
                ## Get convergence error in respect to basis function     
                for N in baseList:
                    N = int(N)
                    if feelpp.Environment.isMasterRank():
                        print("\n\n-----------------------------")
                        print(f"[NIRB] Test with N = {N}, lambda = {lmd}")
                    online_error_sampling(nirb_on, RESPATH,errorfile=errorfile, Nb=N, Nsample=Nsample, Xi_test=Xi_test, verbose=False, regulParam=lmd)
                    del nirb_on.RectificationMat[N] 
        else :
            ## Get convergence error in respect to basis function     
            for N in baseList:
                N = int(N)
                if feelpp.Environment.isMasterRank():
                    print("\n\n-----------------------------")
                    print(f"[NIRB] Test with N = {N}, lambda = {lmd}")
                online_error_sampling(nirb_on, RESPATH,errorfile=errorfile, Nb=N, Nsample=Nsample, Xi_test=Xi_test, verbose=False, regulParam=lmd)
                del nirb_on.RectificationMat[N]
                
    if timeExec:
        ## Get online time mesure
        resOnline = online_time_measure(nirb_on, Nsample=Nsample)
        if feelpp.Environment.isMasterRank():
            res = resOnline
            res['nirb_offline'] = resOffline['nirb_offline']
            res['nproc'] = size
            res['stype'] = idmodel 
            df = pd.DataFrame(res) 
            # save file 
            parent3 = list(Path(RESPATH).absolute().parents)[3]
            file = Path(f"{parent3}/nirb_parallel_time_exec.csv")
            header = not os.path.isfile(file)
            df.to_csv(file, mode='a', index=False, header=header)
            print(f"[NIRB] Execution time mesure are saved in {file}")

    if feelpp.Environment.isMasterRank():
        print("=============================================")
        print(f"Run test_perf_nirb done, doRectification : {doRectification}, doGreedy : {doGreedy}")
        print(f"Result are stored in : {os.path.abspath(RESPATH)}")
        print(f"Set convergence study = {convergence}")
        print(f"Set execution time = {timeExec}")
        print(f"Elapsed time offline = {resOffline}")
        print("=============================================")