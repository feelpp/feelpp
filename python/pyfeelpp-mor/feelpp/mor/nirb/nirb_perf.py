import numpy as np
from tqdm import tqdm
from feelpp.mor.nirb import nirbOffline, nirbOnline
from petsc4py import PETSc


"""
Check convergence
"""
def checkConvergence(nirb_offline: nirbOffline, nirb_online: nirbOnline, Ns=10, Xi_test=None):
    """ Check the convergence of the offline step

    Parameters:
    ----------
        Ns (int) : number of random mu to test
        Xi_test (list) : list of mu to test

    Returns:
    --------
        Error (dict) : dictionnary with the following keys :
            - 'N' : sizes the reduced basis
            - 'l2(uh-uHn)' : l2 norm of the error between the fine and coarse solution
            - 'l2(uh-uHn)rec' : l2 norm of the error between the fine and coarse solution with rectification
            - 'l2(uh-uhn)' : l2 norm of the error between the fine and reduced solution
            - 'l2(uh-uH)' : l2 norm of the error between the u_H^calN interpolated on the fine mesh, and the fine solution

    """
    if nirb_offline.worldcomm.isMasterRank():
        print(f"[NIRB::checkConvergence] Compute offline convergence error")

    if Xi_test is None:
        s = nirb_offline.Dmu.sampling()
        s.sampling(Ns, 'log-random')
        vector_mu = s.getVector()
    else:
        vector_mu = Xi_test

    Ntest = len(vector_mu)

    Nbasis = nirb_offline.N
    assert Nbasis == nirb_online.N, "The number of basis functions in the offline and online phase are not the same"

    fineSolutions= []                   # fine solutions u_h^\N(mu)
    interpolatedSolutions = []          # coarse solutions interpolated in fine mesh u_Hh^\N(mu)
    for j in tqdm(range(len(vector_mu)), desc=f"[NIRB::checkConvergence] Compute fine solutions", ascii=False, ncols=120):
        mu = vector_mu[j]
        uH_interpolated = nirb_online.getInterpolatedSolution(mu)
        fineSolutions.append(nirb_offline.getToolboxSolution(nirb_offline.tbFine, mu))
        interpolatedSolutions.append(uH_interpolated)

    Uhn = nirb_offline.Xh.element()

    Error = {'N':[], 'idx':[], 'l2(uh-uHn)':[], 'l2(uh-uHn)rec':[], 'l2(uh-uhn)' : [], 'l2(uh-uH)':[]}

    pas = 1 if (Nbasis < 50) else 5
    for size in tqdm(range(1, Nbasis + 1, pas), desc=f"[NIRB::checkConvergence] Compute convergence error", ascii=False, ncols=120):

        for j in range(Ntest):

            Uhn         = nirb_online.fineProjection(fineSolutions[j], size) # get reduced sol in a basis function space
            UnirbNoRect = nirb_online.getOnlineSolution(mu, Nb=size, interSol=interpolatedSolutions[j], doRectification=False)
            UnirbRect   = nirb_online.getOnlineSolution(mu, Nb=size, interSol=interpolatedSolutions[j], doRectification=True)

            UnirbNoRect.add(-1, fineSolutions[j])       # UnirbNoRect = u_h^Ncal(mu) - u_Hh^Ncal(mu)
            UnirbRect.add(-1, fineSolutions[j])         # UnirbRect   = u_h^Ncal(mu) - Ru_Hh^Ncal(mu)
            Uhn.add(-1, fineSolutions[j])               # Uhn         = u_h^Ncal(mu) - u_h^\N(mu)
            UhInt = fineSolutions[j] - interpolatedSolutions[j]

            nirbError           = np.sqrt(abs(nirb_offline.l2ScalarProductMatrix.energy(UnirbNoRect, UnirbNoRect)))
            nirbErrorRect       = np.sqrt(abs(nirb_offline.l2ScalarProductMatrix.energy(UnirbRect, UnirbRect)))
            fineProjectionError = np.sqrt(abs(nirb_offline.l2ScalarProductMatrix.energy(Uhn, Uhn)))
            interpolationError  = np.sqrt(abs(nirb_offline.l2ScalarProductMatrix.energy(UhInt, UhInt)))

            Error['l2(uh-uHn)'].append(nirbError)
            Error['l2(uh-uHn)rec'].append(nirbErrorRect)
            Error['l2(uh-uhn)'].append(fineProjectionError)
            Error['l2(uh-uH)'].append(interpolationError)
            Error['N'].append(size)
            Error['idx'].append(j)

    return Error

def computeError(nirb_on: nirbOnline, mu, Nb=-1, h1=False, relative=True):
    """Compute the errors for a given parameter mu

    Parameters
    ----------
    nirb_on : nirbOnline
        initialized nirbOnline class
    mu : parameterSpaceElement
        parameter to test
    Nb : int, optional
        size of the reduced space, by default -1. If -1, the whole basis is used
    h1 : bool, optional
        compute error in h1 norm too, by default False
    relative : bool, optional
        compute relative errors, by default True

    Returns
    -------
    dicrionnary
        see computeErrorSampling
    """
    return computeErrorSampling(nirb_on, Nb=Nb, Xi_test=[mu], h1=h1, relative=relative)

def computeErrorSampling(nirb_on: nirbOnline, Nb=-1, Nsample=1, Xi_test=None, samplingType='log-random', h1=False, relative=True):
    """Compute the convergence errors of the projection of FE solution (in fine mesh) in the reduced space
        and the nirb solution with and without rectification for a given sampling of parameter.
        The rectification matrix has to be defined (so make sure to set doRectification = True in offline and online phase)

    Args:
    -----
        nirb_on (class): initialized nirbOnline class
        Nb (int, optional) : Size of reduced space, by default -1. If -1, the whole basis is used
        Nsample (int, optional) : number of parameter to test, by default 1
        Xi_test (list of ParameterSpaceElement, optional) : list of parameter to test, by default None. If None, a random sampling is done
        samplingType (str, optional) : type of sampling, by default 'log-random'.
        relative (bool) : if True, compute the relative error (error / norm of the true solution) default to True

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
    print(f"[NIRB Perf] ComputeErrorSampling with basis size {Nb} on a sample of {Nsample} parameters")
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

    error = {'idx': [], 'l2(uh-uHn)':[],'l2(uh-uHn)rec':[],'l2(uh-uhn)':[],'l2(uh)':[], 'l2(uh-uH)':[]}
    if h1:
        error.update({'h1(uh-uHn)':[],'h1(uh-uHn)rec':[],'h1(uh-uhn)':[],'h1(uh)':[], 'h1(uh-uH)':[]})

    # for i, mu in enumerate(mus):
    for i, mu in enumerate(tqdm(mus,desc=f"[NIRB] ComputeErrorSampling", ascii=False, ncols=120)):

        uH_interpolated = nirb_on.getInterpolatedSolution(mu)
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)

        uNH  = nirb_on.getOnlineSolution(mu, Nb=Nb, interSol=uH_interpolated, doRectification=False)
        uNHr = nirb_on.getOnlineSolution(mu, Nb=Nb, interSol=uH_interpolated, doRectification=True)
        uhN  = nirb_on.fineProjection(uh, Nb=Nb)

        # error
        error['idx'].append(i)

        rel = nirb_on.normMat(uh, l2Mat) if relative else 1.
        error['l2(uh-uHn)'].append(nirb_on.normMat(uNH - uh, l2Mat) / rel)
        error['l2(uh-uHn)rec'].append(nirb_on.normMat(uNHr - uh, l2Mat) / rel)
        error['l2(uh-uhn)'].append(nirb_on.normMat(uhN - uh, l2Mat) / rel)
        error['l2(uh)'].append(nirb_on.normMat(uh, l2Mat) / rel)
        error['l2(uh-uH)'].append(nirb_on.normMat(uhN - uh,l2Mat) / rel)

        if h1:
            rel = nirb_on.normMat(uh, h1Mat) if relative else 1.
            error['h1(uh-uHn)'].append(nirb_on.normMat(uNH - uh, h1Mat) / rel)
            error['h1(uh-uHn)rec'].append(nirb_on.normMat(uNHr - uh, h1Mat) / rel)
            error['h1(uh-uhn)'].append(nirb_on.normMat(uhN - uh, h1Mat) / rel)
            error['h1(uh)'].append(nirb_on.normMat(uh, h1Mat) / rel)
            error['h1(uh-uH)'].append(nirb_on.normMat(uhN - uh, h1Mat) / rel)

    return error


def computeErrorsH(nirb_on, tbRef, mu, path=None, name=None):
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
