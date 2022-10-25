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

    args = parser.parse_args()
    config_file = args.config_file

    cfg = feelpp.readcfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']
    toolboxType = config_nirb['toolboxType']

    start= time()
    
    nirb_on = nirbOnline(**config_nirb)

    mu = nirb_on.Dmu.element()
    nirb_on.loadData()
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
    errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample)
    
    df = pd.DataFrame(errorN)
    df['N'] = nirb_on.N

    if doRectification:
        file =f"errors{Nsample}ParamsRectif.csv"
    else:
        file =f"errors{Nsample}Params.csv"

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
        



    # ### Solve in very fine mesh 
    # msh_path=f"/data/home/elarif/squareref.msh" # link to a very fine mesh as reference 
    # tbRef = ToolboxModel(dim, H, h, toolboxType, cfg_path, model_path, msh_path, order)

    # mu_ref = tbRef.Dmu.element()
    # mu_ref.setParameters(tbRef.model['Parameters']) # Get defaults parameters from json file 
    
    # """ We can solve the fine problem once and save the solution in a file to load every changement of mesh size 
    #     or recompute this solution in each iteration (costly) """
    # # uh = tbRef.getToolboxSolution(tbRef.tbFine, mu)
    # path = f"nirb_ref_{toolboxType}_sol"
    # name = 'FE_sol_ref'
    # # uh.save(path,'FE_sol_ref')

    # errorH = ComputeErrorsH(nirb_on, tbRef, mu_ref, path, name) # With solution saved in <path>/<name>.h5 
    # # errorH = ComputeErrorsH(nirb_on, tbRef, mu_ref) # With computing the solution 

    # print(f'[NIRB online] mesh size, error = ', errorH)
    

    """ 
    # errors = np.zeros((Nsample, 3))
    errors = []
    for i, mu in enumerate(mus):
        nirb_on.getInterpSol(mu)
        nirb_on.getCompressedSol()
        uHhN = nirb_on.getOnlineSol()
        online1 = nirb_on.onlineSol.to_petsc().vec()[:] # en commentant cette ligne ça produite des nan à la solution onlineSol après computeErrors

        err = nirb_on.computeErrors()
        errors.append(err[:5])

        print(err)
        if not(err[1] > 1e-10 or err[2] > 1e-10 ):
            errors.append(err[:5])
        # errors[i,:] = nirb_on.computeErrors()
        
        # online2 = nirb_on.onlineSol.to_petsc().vec()[:]
    
    errors = np.array(errors)
    print(errors)
    print("\n")
    error_min = errors.min(axis=0)
    error_mean = errors.mean(axis=0)
    error_max = errors.max(axis=0)

    error = np.concatenate((error_min, error_mean[1:], error_max[1:]))
    """


    # print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")