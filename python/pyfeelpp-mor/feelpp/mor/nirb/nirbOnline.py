from inspect import Parameter
from time import time
from feelpp.mor.nirb.nirb import *
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import sys 
import pandas as pd 
from pathlib import Path


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
    error.append((uHh - uh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))
    error.append((uH - uh).l2Norm())
    error.append((uH - uh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))

    return error 

def ComputeErrorsH(nirb_on, tbRef, mu, path=None, name=None):
    """Compute the error between nirb solution and a reference one given from 
           a very fine mesh in respect to mesh size 

    Args:
        nirb_on (class nirbOnline, initialized): nirb online instance 
        tbRef (class ToolBoxModel, initialize): reference toolbox 
        mu (ParameterSpaceElement): parameter 
        path (str) : a path to load the reference solution
        name (str) : the name of the reference solution 
    """

    uHh, _ = nirb_on.getOnlineSol(mu) # nirb solution 
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu) # FE solution 

    if (path != None) and (name != None): 
        uref = tbRef.Xh.element()
        uref.load(path,name)
    else :
        uref = tbRef.getToolboxSolution(tbRef.tbFine, mu) # get reference sol 

    Op = tbRef.createInterpolator(tbRef.tbFine, nirb_on.tbFine)
    Uref_int = Op.interpolate(uref)

    error = [nirb_on.tbFine.mesh().hAverage()]
    error.append((Uref_int - uHh).l2Norm())
    error.append((Uref_int - uHh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))
    error.append((Uref_int - uh).l2Norm())
    error.append((Uref_int - uh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))

    return error

def ComputeErrorSampling(nirb_on, Nsample=1, samplingType='log-random'):
    """Compute the errors between nirb solution and FE one for 
        a given number of parameters 
    Args:
        nirb_on (class): initialized nirbOnline class 
        Nsample (int): number of parameter. Defaults to 1.  
        samplingType (str, optional): type of sampling distribution. Defaults to 'log-random'.

    return : 
    dict: containing
        - dict['mu'] : the index of the parameter from 0... to Nsample
        - dict['l2u-uHn'] : the l2 norm between FE solution (u) and the nirb solution (uHn)
        - dict['lINFu-uHn'] : the infinity norm between FE solution (u) and the nirb solution (uHn)
        - dict['l2u-uH'] : the l2 norm between both FE solution (u in the fine mesh) and (uH in the coarse mesh)
        - dict['lINFu-uH'] : the infinity norm between both FE solution (u in the fine mesh) and (uH in the coarse mesh) 
         
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

    errors = {}
    errors['parameter'] = [*range(Nsample)] 
    errors['l2u-uHn'] = err[:,0]
    errors['lINFu-uHn'] = err[:,1]
    errors['l2u-uH'] = err[:,2]
    errors['lINFu-uH'] = err[:,3]

    return errors


if __name__ == "__main__":

    dim = 3
    order = 1
    if dim == 2:
        H = 0.1  # CoarseMeshSize
        h = H**2 # Fine mesh size
    else:
        # fineness of two grids
        H = 0.5  # CoarseMeshSize
        h = 0.1  # Fine mesh size

    PWD = os.getcwd()
    toolboxType='heat'
    modelfile={'heat':'square/square', 'fluid':'lid-driven-cavity/cfd2d'}
    modelsFolder = f"{PWD}/model/"
    cfg_path = f"{modelsFolder}{modelfile[toolboxType]}.cfg"
    geo_path = f"{modelsFolder}{modelfile[toolboxType]}.geo"
    model_path = f"{modelsFolder}{modelfile[toolboxType]}.json"
    if dim == 3:
        cfg_path = f"{modelsFolder}/thermal-fin-3d/thermal-fin.cfg"
        geo_path = f"{modelsFolder}/thermal-fin-3d/fin.geo"
        model_path = f"{modelsFolder}/thermal-fin-3d/thermal-fin.json"

    e = init_feelpp_environment(toolboxType, cfg_path)

    if len(sys.argv)>=2: 
        H = float(sys.argv[1])
        h = H**2
    

    doRectification=False

    nirb_on = nirbOnline(dim, H, h, toolboxType, cfg_path, model_path, geo_path, order=order, doRectification=doRectification)

    mu = nirb_on.Dmu.element() 
    nirb_on.loadData()
    # uHh = nirb_on.getOnlineSol(mu)


    Nsample = 50
    # error1 = ComputeErrors(nirb_on, mu)
    errorN = ComputeErrorSampling(nirb_on, Nsample=Nsample)

    df = pd.DataFrame(errorN)
    if doRectification:
        file =Path(f"errorParamsRectif/errors{nirb_on.N}.csv")
    else:
        file =Path(f"errorParams/errors{nirb_on.N}.csv")

    file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(file, index=False)

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
    Dmu = nirb_on.Dmu
    s = Dmu.sampling()
    Ns = 1
    s.sampling(Ns, 'log-random')
    mus = s.getVector()

    # err = np.zeros(Nsample)
    # start = time()
    # for i,mu in enumerate(mus):
    #     uH = nirb_on.solveOnline(mu)
    #     uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    #     err[i] = (uH-uh).l2Norm()
    # finish = time()

    time_toolbox_start = time()
    for mu in mus:
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    time_toolbox_finish = time()
    time_toolbox_elapsed = (time_toolbox_finish - time_toolbox_start) / Ns


    time_nirb_start = time()
    for mu in mus:
        uHh = nirb_on.getOnlineSol(mu)
    time_nirb_finish = time()
    time_nirb_elapsed = (time_nirb_finish - time_nirb_start) / Ns

    WriteVecAppend('nirbOnline_time_exec.dat', [nirb_on.N, time_toolbox_elapsed, time_nirb_elapsed])
    print(f"[NIRB online] Elapsed toolbox_time = {time_toolbox_elapsed} nirb_time = {time_nirb_elapsed}")
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
        file1='nirb_error_rectif1.dat'
        fileN=f'nirb_error_rectif{Nsample}.dat'
        fileH=f'nirb_error_rectifH.dat'
    else :
        file1='nirb_error1.dat'
        fileN=f'nirb_error{Nsample}.dat'
        fileH=f'nirb_errorH.dat'
        
    # WriteVecAppend(file1,error1)
    # WriteVecAppend(fileN,errormean)
    # WriteVecAppend(fileH,errorH)

    # print(f"[NIRB] Online Elapsed time =", finish-start)
    print(f"[NIRB] Online part Done !!")