from time import time
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment


def ComputeErrors(nirb_on, mu):
    """Compute the error between nirb solution and FE one for 
        a given parameter 
    Args:
        nirb_on (class): initialized nirbOnline class 
        mu (ParameterSpaceElement) : parameter 

    return : 
    list: containing
        - the number reduced basis N
        - the errors L2 and Linf between nirb solution and FE fine solution
        - the errors L2 and Linf between inteprolated FE coarse solution and FE fine solution
    list : 
    """

    uHh, _ = nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    
    error = [nirb_on.N]

    error.append((uHh - uh).l2Norm())
    error.append((uHh - uh).linftyNorm())
    error.append((uH - uh).l2Norm())
    error.append((uH - uh).linftyNorm())

    return error 

def ComputeErrorSampling(nirb_on, Nsample=1, samplingType='log-random'):
    """Compute the mean error between nirb solution and FE one for 
        a given number of parameter 
    Args:
        nirb_on (class): initialized nirbOnline class 
        Nsample (int): number of parameter. Defaults to 1.  
        samplingType (str, optional): type of sampling distribution. Defaults to 'log-random'.

    return : 
    list: containing
        - the number reduced basis N
        - the errors L2 and Linf between nirb solution and FE fine solution
        - the errors L2 and Linf between inteprolated FE coarse solution and FE fine solution
    list : 
    """

    # finish = timeit()
    Dmu = nirb_on.Dmu
    s = Dmu.sampling()
    s.sampling(Nsample, samplingType)
    mus = s.getVector()

    err = np.zeros((Nsample,4))
    for i,mu in enumerate(mus):
        uH = nirb_on.getInterpSol(mu)
        uHh, _ = nirb_on.getOnlineSol(mu)
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
        
        err[i,0] = (uHh - uh).l2Norm()
        err[i,1] = (uHh - uh).linftyNorm()
        err[i,2] = (uH - uh).l2Norm()
        err[i,3] = (uH - uh).linftyNorm()

    err_mean = err.mean(axis=0)
    error = [nirb_on.N]
    error.append(err_mean[0])
    error.append(err_mean[1])
    error.append(err_mean[2])
    error.append(err_mean[3])

    return error


if __name__ == "__main__":

    # fineness of two grids
    H = 0.1  # CoarseMeshSize
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
    
    # start = time()

    doRectification=False

    nirb_on = nirbOnline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, doRectification=doRectification)

    mu = nirb_on.Dmu.element() 
    nirb_on.loadData()
    # uHh = nirb_on.getOnlineSol(mu)

    Nsample = 10

    error1 = ComputeErrors(nirb_on, mu)
    errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample)
    print("[NIRB online] error 1 param = ", error1)
    print(f"[NIRB online] error {Nsample} param = ", errorN)
    


    """
    err = np.zeros(Nsample)
    start = time()
    for i,mu in enumerate(mus):
        uH = nirb_on.solveOnline(mu)
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
        err[i] = (uH-uh).l2Norm()
    finish = time()

    time_toolbox_start = time()
    for mu in mus:
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
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

    # finish = time() 

    if doRectification:
        file='nirb_error_rectif.dat'
    else :
        file='nirb_error.dat'
    WriteVecAppend(file,error1)

    # print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")