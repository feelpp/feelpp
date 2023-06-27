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
    print(f"[CONVCONV] Check convergence with {Xi_test}")
    if nirb_offline.worldcomm.isMasterRank():
        print(f"[NIRB] Compute offline convergence error")

    if Xi_test is None:
        s = nirb_offline.Dmu.sampling()
        s.sampling(Ns, 'log-random')
        vector_mu = s.getVector()
    else:
        vector_mu = Xi_test

    Ntest = len(vector_mu)

    Nbasis = nirb_offline.N
    assert Nbasis == nirb_online.N, "The number of basis functions in the offline and online phase are not the same"

    print(f"[CONVCONV] Number of basis functions : {Nbasis}")

    interpolationOperator = nirb_offline.createInterpolator(domain=nirb_online.tbCoarse, image=nirb_offline.tbFine)

    fineSolutions= []                   # fine solutions u_h^\N(mu)
    interpolatedSolutions = []          # coarse solutions interpolated in fine mesh u_Hh^\N(mu)
    for mu in vector_mu:
        fineSolutions.append(nirb_offline.getToolboxSolution(nirb_offline.tbFine, mu))
        uH = nirb_offline.getToolboxSolution(nirb_online.tbCoarse, mu)
        interpolatedSolutions.append(interpolationOperator.interpolate(uH))

    Uhn = nirb_offline.Xh.element()

    Error = {'N':[], 'idx':[], 'l2(uh-uHn)':[], 'l2(uh-uHn)rec':[], 'l2(uh-uhn)' : [], 'l2(uh-uH)':[]}

    pas = 1 if (Nbasis < 50) else 5
    for size in tqdm(range(1, Nbasis + 1, pas), desc=f"[NIRB] Compute convergence error", ascii=False, ncols=100):
        print(f"[CONVCONV] Size of reduced basis : {size}")


        for j in range(Ntest):
            Uhn.setZero()

            for k in range(size): # get reduced sol in a basis function space
                projectionCoefficientFine = nirb_offline.l2ScalarProductMatrix.energy(fineSolutions[j], nirb_offline.reducedBasis[k])   # <u_h^Ncal(mu_j), phi_k>L2
                Uhn.add(float(projectionCoefficientFine), nirb_offline.reducedBasis[k])

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
    else:
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
    rel = nirb_on.normMat(uh, l2Mat) if relative else 1
    error['l2(uh-uHn)'].append(nirb_on.normMat(uNH-uh, l2Mat) / rel)
    error['l2(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, l2Mat) / rel)
    error['l2(uh-uhn)'].append(nirb_on.normMat(uNh-uh, l2Mat) / rel)
    error['l2(uh)'].append(nirb_on.normMat(uh, l2Mat))
    error['l2(uh-uH)'].append(nirb_on.normMat(uH-uh, l2Mat) / rel)
    if h1:
        rel = nirb_on.normMat(uh, h1Mat) if relative else 1
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
        uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)

        uNH = getNirbProjection(nirb_on, uH, Nb=Nb)
        uNHr = getNirbProjection(nirb_on, uH, doRectification=True, Nb=Nb, regulParam=regulParam)
        uNh = getNirbProjection(nirb_on, uh, Nb=Nb)

        # error
        error['l2(uh-uHn)'].append(nirb_on.normMat(uNH-uh, l2Mat))
        error['l2(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, l2Mat))
        error['l2(uh-uhn)'].append(nirb_on.normMat(uNh-uh, l2Mat))
        error['l2(uh)'].append(nirb_on.normMat(uh, l2Mat))
        error['l2(uh-uH)'].append(nirb_on.normMat(uH-uh,l2Mat))
        if h1:
            error['h1(uh-uHn)'].append(nirb_on.normMat(uNH-uh, h1Mat))
            error['h1(uh-uHn)rec'].append(nirb_on.normMat(uNHr-uh, h1Mat))
            error['h1(uh-uhn)'].append(nirb_on.normMat(uNh-uh, h1Mat))
            error['h1(uh)'].append(nirb_on.normMat(uh, h1Mat))
            error['h1(uh-uH)'].append(nirb_on.normMat(uH-uh, h1Mat))

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

def getNirbProjection(nirb_on, u, Nb=None, doRectification=False, regulParam=1.-10):
    """Get the projection of a given discrete function in the reduced space with or without rectification

    Args
    ----
        nirb_on (class): initialized nirbOnline class
        u (feelpp._discr.element) : function
        Nb (int, optional) : Size of reduced space, by default None. If None, the whole basis is used
        doRectification (bool) : default to True
        regulParam(float, optionnal) : regularization parameter of the rectification pre-process. Defaults to 1.e-10
    """
    if Nb is None : Nb=nirb_on.N
    if doRectification:
        assert nirb_on.doRectification, f"set doRectification from nirb_on class"

    coef = nirb_on.getCompressedSol(solution=u, Nb=Nb)

    uNh = nirb_on.Xh.element()
    uNh.setZero()

    if doRectification:
        if Nb not in nirb_on.RectificationMat :
                nirb_on.RectificationMat[Nb] = nirb_on.getRectification(nirb_on.coeffCoarse, nirb_on.coeffFine, Nb=Nb, lambd=regulParam)
        coef = nirb_on.RectificationMat[Nb] @ coef

    for i in range(Nb):
        uNh.add(float(coef[i]), nirb_on.reducedBasis[i])

    return uNh
