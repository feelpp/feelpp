from distutils.log import error
from timeit import timeit
from time import time
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment

if __name__ == "__main__":

    # fineness of two grids
    H = 0.001  # CoarseMeshSize
    h = H**2 # fine mesh size
    dim = 2

    PWD = os.getcwd()
    toolboxType='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxType]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxType]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxType]}.json"

    e = init_feelpp_environment(toolboxType, cfg_path)
    
    start = time()

    doRectification=False

    nirb_on = nirbOnline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)

    nirb_on.loadData()

    ## Compute NIRB error with one parameter 
    """ 
    # nirb_on.getInterpSol()
    # nirb_on.getCompressedSol()
    # nirb_on.getOnlineSol()

    # online1 = nirb_on.onlineSol.to_petsc().vec()[:] # en commentant cette ligne ça produite des nan à la solution onlineSol après computeErrors
  
    # error = nirb_on.computeErrors()
    
    # online2 = nirb_on.onlineSol.to_petsc().vec()[:]

    # print("[nirb main] diff online sol befor/after error", np.max(np.abs(online1-online2)))

    # print("----------- Errors ----------------------")
    # print(f"Nb mode = {error[0]}")
    # print(f"L2 norm= {error[1]}")
    # print(f"Inf norm = {error[2]}")
    # print('error =', error)
    """


    # finish = timeit()
    Dmu = nirb_on.Dmu
    Nsample = 1
    s = Dmu.sampling()
    s.sampling(Nsample, 'log-random')
    mus = s.getVector()
    ## Compute mean error of a sampling of parameters 
    err = np.zeros((Nsample,4))
    for i,mu in enumerate(mus):
        uH = nirb_on.getInterpSol(mu)
        nirb_on.getCompressedSol()
        uHh = nirb_on.getOnlineSol()
        print(i, 'ok')
        online1 = nirb_on.onlineSol.to_petsc().vec()[:]
        uh = nirb_on.getSolution(nirb_on.tbFine, mu)
        err[i,0] = (uHh - uh).l2Norm()
        err[i,1] = (uHh - uh).linftyNorm()
        err[i,2] = (uH - uh).l2Norm()
        err[i,3] = (uH - uh).linftyNorm()
        print(i, 'okk')
    err_mean = err.mean(axis=0)
    error = [nirb_on.N]
    error.append(err_mean[0])
    error.append(err_mean[1])
    error.append(err_mean[2])
    error.append(err_mean[3])
    
    print("[NIRB online] error = ", error)

        
        
    """
    err = np.zeros(Nsample)
    start = time()
    for i,mu in enumerate(mus):
        uH = nirb_on.solveOnline(mu)
        uh = nirb_on.getSolution(nirb_on.tbFine, mu)
        err[i] = (uH-uh).l2Norm()
    finish = time()

    time_toolbox_start = time()
    for mu in mus:
        uh = nirb_on.getSolution(nirb_on.tbFine, mu)
    time_toolbox_finish = time()
    time_toolbox_elapsed = (time_toolbox_finish - time_toolbox_start) / Nsample


    time_nirb_start = time()
    for mu in mus:
        nirb_on.getInterpSol(mu)
        nirb_on.getCompressedSol()
        uHh = nirb_on.getOnlineSol()
    time_nirb_finish = time()
    time_nirb_elapsed = (time_nirb_finish - time_nirb_start) / Nsample

    WriteVecAppend('nirbOnline_time_exec.dat', [nirb_on.N, time_toolbox_elapsed, time_nirb_elapsed])
    """
   
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

    finish = time() 

    if doRectification:
        file='nirb_error_rectif.dat'
    else :
        file='nirb_error.dat'
    WriteVecAppend(file,error)

    print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")