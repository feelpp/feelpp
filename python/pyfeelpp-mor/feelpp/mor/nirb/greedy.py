# -*- coding: utf-8 -*-
## functions for greedy algorithm in NIRB
## Thomas Saigre, Ali Elarif
## 01/2023

from tqdm import tqdm
from feelpp.mor.nirb import nirbOffline, nirbOnline

def initProblemGreedy(offline: nirbOffline, online: nirbOnline, Ninit, Ntrain, tol=1e-5, Xi_train=None, Nmax=50, samplingMode="log-random", Mexpr='N-1'):
    """Generate the reduced space using the greedy algorithm.

    Parameters
    ----------
    offline : nirbOffline
        nirbOffline object, already initialized
    online : nirbOnline
        nirbOnline object, already initialized
    Ninit : int
        size of the initial reduced space
    Ntrain : int
        size of the train parameter set
    tol : float, optional
        tolerance for greedy algorithm, by default 1e-5
    Xi_train : list, optional
        train set for the greedy algorithm, by default None. If None, a sampling of size Ntrain is generated
    Nmax : int, optional
        maximal size of the basis, by default 50
    samplingMode : str, optional
        sampling mode for Xi_train, by default "log-random"
    Mexpr : str, optional
        expression of how M is computed regarding N, by default 'N-1'

    Returns
    -------
    tuple
        TODO
    """

    if offline.tbCoarse is None:
        raise Exception("Coarse toolbox needed for computing coarse Snapshot. set initCoarse->True in initialization")
    Nmax = min(Nmax, Ntrain)

    if Xi_train is None:
        s = offline.Dmu.sampling()
        s.sampling(20, samplingMode)
        Xi_train = s.getVector()
    Xi_train_copy = Xi_train.copy()

    Delta_star = tol + 1
    Deltas_conv = []
    S = []
    offline.fineSnapShotList = []
    offline.coarseSnapShotList = []

    # Computation of coarse solutions
    coarseSolutions = {}
    interpolatedSolutions = {}
    for mu in tqdm(Xi_train_copy, desc="[NIRB] Computing interpolated coarse solutions", ncols=120):
        coarseSolutions[mu] = offline.getToolboxSolution(offline.tbCoarse, mu)
        interpolatedSolutions[mu] = online.getInterpolatedSolution(mu)

    # intialization of the first reduced space
    for i in range(Ninit):
        mu0 = offline.Dmu.element()
        S.append(mu0)
        offline.fineSnapShotList.append(   offline.getToolboxSolution(offline.tbFine, mu0)    )
        offline.coarseSnapShotList.append( offline.getToolboxSolution( offline.tbCoarse, mu0) )

    RIC = offline.generateReducedBasis()
    print("[NIRB] RIC = ", RIC)
    online.setBasis(offline)
    N = Ninit
    assert N == offline.N

    # Greedy loop
    if offline.worldcomm.isMasterRank():
        print(f"[NIRB] Starting greedy loop with tolerance tol = {tol} and Nmax = {Nmax}")
    while Delta_star > tol and N < Nmax:
        if offline.worldcomm.isMasterRank():
            print(f"[NIRB] Greedy loop : N = {N} Delta_star = {Delta_star}")
        M = eval(Mexpr)

        mu_star, Xi_train, Delta_star = greedyStep(offline, online, Xi_train, N, M, interpolatedSolutions)
        N = offline.N
        Deltas_conv.append(Delta_star)
        S.append(mu_star)


    online.saveExporter()

    if offline.worldcomm.isMasterRank():
        print(f"[NIRB] Number of snapshot computed : {N}")
        print(f"[NIRB] Convergence of the greedy algorithm : {Deltas_conv}")

    return S, Xi_train


def greedyStep(offline: nirbOffline, online: nirbOnline, Xi_train, N, M, interpSol):
    """Compute on step of greedy algorithm

    Parameters
    ----------
    offline : nirbOffline
        nirbOffline object
    online : nirbOnline
        nirbOnline object
    Xi_train : list of ParameterSpaceElement
        train set of parameters
    N : int
        size of reduced basis
    M : int
        size of rediced basis to compare
    interpSol : _type_
        Interpolated solution

    Returns
    -------
    tuple
        TODO
    """

    Delta_star = -float('inf')

    Deltas = []

    for i, mu in enumerate(tqdm(Xi_train,desc=f"[NIRB] Greedy selection", ascii=False, ncols=120)):
        uHhN = online.getOnlineSolution(mu, N, interSol=interpSol[mu])
        uHhM = online.getOnlineSolution(mu, M, interSol=interpSol[mu])

        diff = uHhN - uHhM
        Delta = offline.l2ScalarProductMatrix.energy(diff, diff)

        Deltas.append(Delta)

        if Delta > Delta_star:
            Delta_star = Delta
            mu_star = mu
            idx = i
    uHN_export = online.getOnlineSolution(mu_star, N, interSol=interpSol[mu_star])
    uHM_export = online.getOnlineSolution(mu_star, M, interSol=interpSol[mu_star])
    uh = online.getToolboxSolution(online.tbFine, mu_star)

    if online.exporter is not None:
        online.exportField(uHN_export, f"NIRB{N:02}uHhN")
        online.exportField(uHM_export, f"NIRB{N:02}uHhM")
        online.exportField(uh, f"NIRB{N:02}uh")

    del Xi_train[idx]
    uH = offline.getToolboxSolution( offline.tbCoarse, mu_star)
    uh = offline.getToolboxSolution( offline.tbFine, mu_star)
    offline.addFunctionToBasis( uh, coarseSnapshot=uH )
    online.setBasis(offline)

    return mu_star, Xi_train, Delta_star
