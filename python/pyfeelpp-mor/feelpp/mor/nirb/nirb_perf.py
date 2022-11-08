import numpy as np
from petsc4py import PETSc


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
    l2Mat, h1Mat = nirb_on.generateOperators(h1=True)

    uHh= nirb_on.getOnlineSol(mu)
    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    
    error = [nirb_on.N]

    error.append(nirb_on.normMat(l2Mat, uHh-uh)) # L2 norm 
    error.append(nirb_on.normMat(h1Mat, uHh-uh)) # H1 norm 
    error.append(nirb_on.normMat(l2Mat, uH-uh))
    error.append(nirb_on.normMat(h1Mat, uH-uh))

    # error.append((uHh - uh).l2Norm())
    # error.append((uHh - uh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))
    # error.append((uH - uh).l2Norm())
    # error.append((uH - uh).to_petsc().vec().norm(PETSc.NormType.NORM_INFINITY))

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

    uHh= nirb_on.getOnlineSol(mu) # nirb solution 
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

def ComputeErrorSampling(nirb_on, Nsample=1, samplingType='log-random', Xi_test=None):
    """Compute the errors between nirb solution and FE one for 
        a given number of parameters 
    Args:
        nirb_on (class): initialized nirbOnline class 
        Nsample (int): number of parameter. Defaults to 1.  
        samplingType (str, optional): type of sampling distribution. Defaults to 'log-random'.

    return : 
    dict: containing
        - dict['parameter'] : the index of the parameter from 0... to Nsample
        - dict['l2(u-uHn)'] : the l2 norm between FE solution (u) and the nirb solution (uHn)
        - dict['h1(u-uHn)'] : the h1 norm between FE solution (u) and the nirb solution (uHn)
        - dict['l2(u-uH)'] : the l2 norm between both FE solution (u in the fine mesh) and (uH in the coarse mesh)
        - dict['h1(u-uH)'] : the h1 norm between both FE solution (u in the fine mesh) and (uH in the coarse mesh)
    """

    # finish = timeit()
    if Xi_test is None:
        Dmu = nirb_on.Dmu
        s = Dmu.sampling()
        s.sampling(Nsample, samplingType)
        mus = s.getVector()
    else:
        Nsample = len(Xi_test)
        mus = Xi_test

    l2Mat, h1Mat = nirb_on.generateOperators(h1=True)

    err = np.zeros((Nsample,4))
    for i,mu in enumerate(mus):
        uH = nirb_on.getInterpSol(mu)
        uHh= nirb_on.getOnlineSol(mu)
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)

        err[i,0] = nirb_on.normMat(l2Mat, uHh-uh)
        err[i,1] = nirb_on.normMat(h1Mat, uHh-uh)
        err[i,2] = nirb_on.normMat(l2Mat, uH-uh)
        err[i,3] = nirb_on.normMat(h1Mat, uH-uh)

    errors = {}
    errors['parameter'] = [*range(Nsample)] 
    errors['l2(u-uHn)'] = err[:,0]
    errors['h1(u-uHn)'] = err[:,1]
    errors['l2(u-uH)'] = err[:,2]
    errors['h1(u-uH)'] = err[:,3]

    return errors