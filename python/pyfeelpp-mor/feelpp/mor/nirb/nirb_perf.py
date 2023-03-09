import numpy as np
from tqdm import tqdm
from petsc4py import PETSc



def ComputeErrors(nirb_on,mu, Nb=None,h1=False, relative=True):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space
        and the nirb solution with and without rectification for a given one parameter.
        The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)

    Args:
    -----
        nirb_on (class): initialized nirbOnline class
        mu (ParameterSpaceElement) : parameter
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used
        h1 (bool) : if True, also compute the h1 norm. Default to False
        relative (bool) : if True, compute the relative error (error / norm of the true solution) default to True

    Returns:
    --------
    dict: containing
        - error['l2(uh-uHn)'] : the L2 norm between True sol and nirb sol without rectification
        - error['l2(uh-uHn)rec'] : the L2 norm between True Sol and nirb sol with rectification
        - error['l2(uh-uhn)'] : the L2 norm between True sol and projection of this True sol in reduced space.
        - error['l2(uh-uH)'] : the L2 norm between True sol and interpolate sol from coarse to fine mesh.
        - error['l2(uh)'] : the L2 norm of the True sol
        N.B : True sol = FE solution in the fine mesh.
              If relativ is True, the error is divided by the norm of the True sol.
              If h1 is True, the same norm are computed in h1 norm.
    """
    if h1: 
        l2Mat, h1Mat = nirb_on.generateOperators(h1=True)
    else :
        l2Mat = nirb_on.generateOperators()

    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    
    uNH = getNirbProjection(nirb_on, uH, Nb=Nb)
    uNHr = getNirbProjection(nirb_on, uH, doRectification=True, Nb=Nb)
    uNh = getNirbProjection(nirb_on, uh, Nb=Nb)

    
    error = {'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})
        
    # error 
    rel = nirb_on.normMat(l2Mat,uh) if relative else 1
    error['l2(uh-uHn)'].append(nirb_on.normMat(uNH-uh, l2Mat) / rel)
    error['l2(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, l2Mat) / rel)
    error['l2(uh-uhn)'].append(nirb_on.normMat(uNh-uh, l2Mat) / rel)
    error['l2(uh)'].append(nirb_on.normMat(uh, l2Mat))
    error['l2(uh-uH)'].append(nirb_on.normMat(uH-uh, l2Mat) / rel)
    if h1 : 
        rel = nirb_on.normMat(h1Mat,uh) if relative else 1
        error['h1(uh-uHn)'].append(nirb_on.normMat(uNH-uh, h1Mat) / rel)
        error['h1(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, h1Mat) / rel)
        error['h1(uh-uhn)'].append(nirb_on.normMat(uNh-uh, h1Mat) / rel)
        error['h1(uh)'].append(nirb_on.normMat(uh, h1Mat))
        error['h1(uh-uH)'].append(nirb_on.normMat(uH-uh, h1Mat) / rel) 


    return error

def ComputeErrorSampling(nirb_on, Nb=None, Nsample = 1, Xi_test=None, samplingType='log-random', h1=False, regulParam=1.e-10):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space
        and the nirb solution with and without rectification for a given sampling of parameter.
        The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)

    Args:
    -----
        nirb_on (class): initialized nirbOnline class
        mu (ParameterSpaceElement) : parameter
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used

    Returns:
    --------
    error (dict) : containing
        - error['l2(uh-uHn)'] : the L2 norm between True sol and nirb sol without rectification
        - error['l2(uh-uHn)rec'] : the L2 norm between True Sol and nirb sol with rectification
        - error['l2(uh-uhn)'] : the L2 norm between True sol and projection of this True sol in reduced space.
        - error['l2(uh-uH)'] : the L2 norm between True sol and interpolate sol from coarse to fine mesh.

        The same keys are given in case of h1 norm with 'h1' instead of 'l2' (exp : 'h1(uh - uHn)' )
        N.B : True sol = FE solution in the fine mesh
    """
    if h1:
        l2Mat, h1Mat = nirb_on.generateOperators(h1=True)
    else :
        l2Mat = nirb_on.generateOperators()

    if Xi_test is None :
        Dmu = nirb_on.Dmu
        s = Dmu.sampling()
        s.sampling(Nsample, samplingType)
        mus = s.getVector()
    else :
        mus = Xi_test

    error = {'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})

    # for i, mu in enumerate(mus): 
    for i, mu in enumerate(tqdm(mus,desc=f"[NIRB] ComputeErrorSampling", ascii=False, ncols=120)):

        uH = nirb_on.getInterpSol(mu)
        # uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
        if not nirb_on.time_dependent:
            uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
        else :
            uh = nirb_on.getTimeToolboxSolution(nirb_on.tbFine, mu)[nirb_on.Ntime-1]
            
        uNH = getNirbProjection(nirb_on, uH, Nb=Nb)
        uNHr = getNirbProjection(nirb_on, uH, doRectification=True, Nb=Nb, regulParam=regulParam)
        uNh = getNirbProjection(nirb_on, uh, Nb=Nb)

        # error 
        error['l2(uh-uHn)'].append(nirb_on.normMat(uNH-uh, l2Mat))
        error['l2(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, l2Mat))
        error['l2(uh-uhn)'].append(nirb_on.normMat(uNh-uh, l2Mat))
        error['l2(uh)'].append(nirb_on.normMat(uh, l2Mat))
        error['l2(uh-uH)'].append(nirb_on.normMat(uH-uh,l2Mat))
        if h1 : 
            error['h1(uh-uHn)'].append(nirb_on.normMat(uNH-uh, h1Mat))
            error['h1(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, h1Mat))
            error['h1(uh-uhn)'].append(nirb_on.normMat(uNh-uh, h1Mat))
            error['h1(uh)'].append(nirb_on.normMat(uh, h1Mat))
            error['h1(uh-uH)'].append(nirb_on.normMat(uH-uh, h1Mat))    

    return error 


def ComputeErrorSamplingTime(nirb_on, Nb=None, Nsample = 1, Xi_test=None, samplingType='log-random', h1=False, regulParam=1.e-10):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space
        and the nirb solution with and without rectification for a given sampling of parameter in a parabolic system. 
        The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)
        The error is computed in the maximum norm in time. 

    Args:
    -----
        nirb_on (class): initialized nirbOnline class
        mu (ParameterSpaceElement) : parameter
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used

    Returns:
    --------
    error (dict) : containing
        - error['l2(uh-uHn)'] : the L2 norm between True sol and nirb sol without rectification
        - error['l2(uh-uHn)rec'] : the L2 norm between True Sol and nirb sol with rectification
        - error['l2(uh-uhn)'] : the L2 norm between True sol and projection of this True sol in reduced space.
        - error['l2(uh-uH)'] : the L2 norm between True sol and interpolate sol from coarse to fine mesh.

        The same keys are given in case of h1 norm with 'h1' instead of 'l2' (exp : 'h1(uh - uHn)' )
        N.B : True sol = FE solution in the fine mesh
    """
    assert nirb_on.time_dependent, f"Use this function only for parabolic problem. Otherwise use ComputeErrorSampling()"

    if h1:
        l2Mat, h1Mat = nirb_on.generateOperators(h1=True)
    else :
        l2Mat = nirb_on.generateOperators()

    if Xi_test is None :
        Dmu = nirb_on.Dmu
        s = Dmu.sampling()
        s.sampling(Nsample, samplingType)
        mus = s.getVector()
    else :
        mus = Xi_test

    error = {'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})

    for i, mu in enumerate(tqdm(mus,desc=f"[NIRB] ComputeErrorSampling", ascii=False, ncols=120)):

        list_uH = nirb_on.getInterpSol(mu,allTimeSol=True)
        list_uh = nirb_on.getTimeToolboxSolution(nirb_on.tbFine, mu)
        
        l2err= np.zeros(5)
        h1err= np.zeros(5)

        for i in range(len(list_uH)):
            uH = list_uH[i]
            uh = list_uh[i]

            uNH = getNirbProjection(nirb_on, uH, Nb=Nb,itr=i)
            uNHr = getNirbProjection(nirb_on, uH, doRectification=True, Nb=Nb, regulParam=regulParam, itr=i)
            uNh = getNirbProjection(nirb_on, uh, Nb=Nb, itr=i)

            l2err[0] = max(l2err[0], nirb_on.normMat(uNH-uh, l2Mat))
            l2err[1] = max(l2err[1], nirb_on.normMat(uNHr-uh, l2Mat))
            l2err[2] = max(l2err[2], nirb_on.normMat(uNh-uh, l2Mat))
            l2err[3] = max(l2err[3], nirb_on.normMat(uh, l2Mat))
            l2err[4] = max(l2err[4], nirb_on.normMat(uH-uh,l2Mat))
            if h1 : 
                h1err[0] = max(h1err[0], nirb_on.normMat(uNH-uh, h1Mat))
                h1err[1] = max(h1err[1], nirb_on.normMat(uNHr-uh, h1Mat))
                h1err[2] = max(h1err[2], nirb_on.normMat(uNh-uh, h1Mat))
                h1err[3] = max(h1err[3], nirb_on.normMat(uh, h1Mat))
                h1err[4] = max(h1err[4], nirb_on.normMat(uH-uh,h1Mat))

        # error 
        error['l2(uh-uHn)'].append(l2err[0])
        error['l2(uh-uHn)rec'].append(l2err[1])
        error['l2(uh-uhn)'].append(l2err[2])
        error['l2(uh)'].append(l2err[3])
        error['l2(uh-uH)'].append(l2err[4])
        if h1 : 
            error['h1(uh-uHn)'].append(h1err[0])
            error['h1(uh-uHn)rec'].append(h1err[1])
            error['h1(uh-uhn)'].append(h1err[2])
            error['h1(uh)'].append(h1err[3])
            error['h1(uh-uH)'].append(h1err[4])    

    return error 


def ComputeErrorsH(nirb_on, tbRef, mu, path=None, name=None):
    """Compute the error between nirb solution and a reference one given from
           a very fine mesh in respect to mesh size

    Args:
    -----
        nirb_on (class nirbOnline, initialized): nirb online instance
        tbRef (class ToolBoxModel, initialize): reference toolbox
        mu (ParameterSpaceElement): parameter
        path (str) : a path to load the reference solution
        name (str) : the name of the reference solution

    Returns:
    --------
    list: containing
        - the number reduced basis N
        - the errors L2 and Linf between nirb solution and FE fine solution
        - the errors L2 and Linf between inteprolated FE coarse solution and FE fine solution
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

def getNirbProjection(nirb_on, u, Nb=None, doRectification=False, regulParam=1.-10, itr=None):
    """Get the projection of a given discrete function in the reduced space with or without rectification

    Args
    ----
        nirb_on (class): initialized nirbOnline class
        u (feelpp._discr.element) : function
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used
        doRectification (bool) : default to True
        regulParam(float, optionnal) : regularization parameter of the rectification pre-process. Defaults to 1.e-10
        itr (int, optional) : the time iteration. Defaults to None. 
    """
    if Nb is None : Nb=nirb_on.Nmu
    if (itr is None) and (nirb_on.time_dependent) : itr = nirb_on.Ntime-1  
    if doRectification:
        assert nirb_on.doRectification, f"set doRectification from nirb_on class"

    coef = nirb_on.getCompressedSol(solution=u, Nb=Nb, itr=itr)

    uNh = nirb_on.Xh.element()
    uNh.setZero()

    if doRectification:
        if Nb not in nirb_on.RectificationMat :
                nirb_on.RectificationMat[Nb] = nirb_on.getRectification(nirb_on.coeffCoarse, nirb_on.coeffFine, Nb=Nb, lambd=regulParam, itr=itr)
        coef = nirb_on.RectificationMat[Nb] @ coef

    for i in range(Nb):
        uNh.add(float(coef[i]), nirb_on.reducedBasis[i])
    
    return uNh 