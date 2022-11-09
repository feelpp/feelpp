import numpy as np
from petsc4py import PETSc



def ComputeErrors(nirb_on, mu, h1=False, relative=True):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space 
        and the nirb solution with and without rectification for a given one parameter. 
         The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)

    Args:
        nirb_on (class): initialized nirbOnline class 
        mu (ParameterSpaceElement) : parameter 

    return : 
    dict: containing
        - dict['N'] : the number of basis function 
        - dict['nirb'] : the L2 norm between True sol and nirb sol without rectification 
        - dict['nirb_r'] : the L2 norm between True Sol and nirb sol with rectification 
        - dict['fine'] : the L2 norm between True sol and projection of this True sol in reduced space. 
        N.B : True sol = FE solution in the fine mesh  
    """
    if h1: 
        l2Mat, h1Mat = nirb_on.generateOperators(h1=True)
    else :
        l2Mat = nirb_on.generateOperators()

    uH = nirb_on.getInterpSol(mu)
    uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
    
    uNH = getNirbProjection(nirb_on, uH)
    uNHr = getNirbProjection(nirb_on, uH, doRectification=True)
    uNh = getNirbProjection(nirb_on, uh)

    
    error = {'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})
        
    # error 
    error['l2(uh-uHn)'].append(nirb_on.normMat(l2Mat, uNH-uh))
    error['l2(uh-uHn)rec'].append(nirb_on.normMat(l2Mat, uNHr-uh))
    error['l2(uh-uhn)'].append(nirb_on.normMat(l2Mat, uNh-uh))
    error['l2(uh)'].append(nirb_on.normMat(l2Mat,uh))
    error['l2(uh-uH)'].append(nirb_on.normMat(l2Mat, uH-uh))
    if h1 : 
        error['h1(uh-uHn)'].append(nirb_on.normMat(h1Mat, uNH-uh))
        error['h1(uh-uHn)rec'].append(nirb_on.normMat(h1Mat, uNHr-uh))
        error['h1(uh-uhn)'].append(nirb_on.normMat(h1Mat, uNh-uh))
        error['h1(uh)'].append(nirb_on.normMat(h1Mat,uh))
        error['h1(uh-uH)'].append(nirb_on.normMat(h1Mat, uH-uh)) 


    return error 

def ComputeErrorSampling(nirb_on, Nsample = 1, samplingType='log-random', h1=False):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space 
        and the nirb solution with and without rectification for a given sampling of parameter. 
         The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)

    Args:
        nirb_on (class): initialized nirbOnline class 
        mu (ParameterSpaceElement) : parameter 

    return : 
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

    Dmu = nirb_on.Dmu
    s = Dmu.sampling()
    s.sampling(Nsample, samplingType)
    mus = s.getVector()

    error = {'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})
        

    for i, mu in enumerate(mus): 

        uH = nirb_on.getInterpSol(mu)
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
        
        uNH = getNirbProjection(nirb_on, uH)
        uNHr = getNirbProjection(nirb_on, uH, doRectification=True)
        uNh = getNirbProjection(nirb_on, uh)

        # error 
        error['l2(uh-uHn)'].append(nirb_on.normMat(l2Mat, uNH-uh))
        error['l2(uh-uHn)rec'].append(nirb_on.normMat(l2Mat, uNHr-uh))
        error['l2(uh-uhn)'].append(nirb_on.normMat(l2Mat, uNh-uh))
        error['l2(uh)'].append(nirb_on.normMat(l2Mat,uh))
        error['l2(uh-uH)'].append(nirb_on.normMat(l2Mat, uH-uh))
        if h1 : 
            error['h1(uh-uHn)'].append(nirb_on.normMat(h1Mat, uNH-uh))
            error['h1(uh-uHn)rec'].append(nirb_on.normMat(h1Mat, uNHr-uh))
            error['h1(uh-uhn)'].append(nirb_on.normMat(h1Mat, uNh-uh))
            error['h1(uh)'].append(nirb_on.normMat(h1Mat,uh))
            error['h1(uh-uH)'].append(nirb_on.normMat(h1Mat, uH-uh))    

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

def ComputeErrorSamplingOld(nirb_on, Nsample=1, samplingType='log-random'):
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
    Dmu = nirb_on.Dmu
    s = Dmu.sampling()
    s.sampling(Nsample, samplingType)
    mus = s.getVector()

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


def getNirbProjection(nirb_on,u,doRectification=False):
    """Get the projection of a given discrete function in the reduced space with or without rectification

    Args:
        nirb_on (class): initialized nirbOnline class
        u (feelpp._discr.element) : function 
        doRectification (bool) : default to True  
    """
    if doRectification:
        assert nirb_on.doRectification, f"set doRectification from nirb_on class"

    coef = nirb_on.getCompressedSol(solution=u)

    uNh = nirb_on.Xh.element()
    uNh.setZero()

    if doRectification:
        coef = nirb_on.RectificationMat@coef
        for i in range(nirb_on.N):
            uNh = uNh + float(coef[i])*nirb_on.reducedBasis[i]
    else :
        for i in range(nirb_on.N):
            uNh = uNh + float(coef[i])*nirb_on.reducedBasis[i]

    return uNh 