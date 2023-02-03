# -*- coding: utf-8 -*-
## utilis function for NIRB
## Thomas Saigre, Ali Elarif
## 01/2023

from tqdm import tqdm


def initProblemGreedy(offline, online, Ninit, Ntrain, eps=1e-5, Xi_train=None, Nmax=50, samplingMode="log-random", Mexpr='N-1'):
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
    eps : float, optional
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
        TO BE COMPLETED
    """

    if offline.tbCoarse is None:
        raise Exception("Coarse toolbox needed for computing coarse Snapshot. set initCoarse->True in initialization")
    Nmax = min(Nmax, Ntrain)

    if Xi_train is None:
        s = offline.Dmu.sampling()
        s.sampling(Ntrain, samplingMode)
        Xi_train = s.getVector()
    Xi_train_copy = Xi_train.copy()

    Delta_star = eps + 1
    Deltas_conv = []
    S = []
    offline.fineSnapShotList = []
    offline.coarseSnapShotList = []

    # Computation of coarse solutions
    coarseSolutions = {}
    interpSol = {}
    for mu in tqdm(Xi_train_copy, desc="[NIRB] Computing interpolated coarse solutions", ncols=120):
        coarseSolutions[mu] = offline.getToolboxSolution(offline.tbCoarse, mu)
        interpSol[mu] = online.getInterpSol(mu)

    # intialization of the first reduced space
    for i in range(Ninit):
        mu0 = offline.Dmu.element()
        S.append(mu0)
        offline.fineSnapShotList.append( offline.getToolboxSolution(offline.tbFine  , mu0) )
        uH = offline.getToolboxSolution( offline.tbCoarse, mu0)
        offline.coarseSnapShotList.append(uH)


    RIC = offline.generateReducedBasis()
    print("[NIRB] RIC = ", RIC)
    online.setBasis(offline)
    N = Ninit
    assert N == offline.N

    # Greedy loop
    if offline.worldcomm.isMasterRank():
        print("[NIRB] Starting greedy loop with tolerance eps = ", eps, " and Nmax = ", Nmax)
    while Delta_star > eps and N < Nmax:
        M = eval(Mexpr)

        mu_star, Xi_train = greedyStep(offline, online, Xi_train, N, M, interpSol)

        S.append(mu_star)      


    online.saveExporter()

    if offline.worldcomm.isMasterRank():
        print(f"[NIRB] Number of snapshot computed : {N}")

    return S, Xi_train


def greedyStep(offline, online, Xi_train, N, M, interpSol):

    Delta_star = -float('inf')

    Deltas = []

    for i, mu in enumerate(tqdm(Xi_train,desc=f"[NIRB] Greedy selection", ascii=False, ncols=120)):
        uHhN = online.getOnlineSol(mu, N, interSol=interpSol[mu])
        uHhM = online.getOnlineSol(mu, M, interSol=interpSol[mu])

        diff = uHhN - uHhM
        Delta = offline.l2ScalarProductMatrix.energy(diff, diff)

        Deltas.append(Delta)

        if Delta > Delta_star:
            Delta_star = Delta
            mu_star = mu
            idx = i

    uHN_export = online.getOnlineSol(mu_star, N, interSol=interpSol[mu_star])
    uHM_export = online.getOnlineSol(mu_star, M, interSol=interpSol[mu_star])
    uh = online.getToolboxSolution(online.tbFine, mu)

    online.exportField(uHN_export, f"NIRB{N:02}uHhN")
    online.exportField(uHM_export, f"NIRB{N:02}uHhM")
    online.exportField(uh, f"NIRB{N:02}uh")

    del Xi_train[idx]
    uH = offline.getToolboxSolution( offline.tbCoarse, mu_star)
    offline.addFunctionToBasis( offline.getToolboxSolution(offline.tbFine, mu_star), coarseSnapshot=uH )
    N += 1
    online.setBasis(offline)

    return mu_star, Xi_train

    if offline.worldcomm.isMasterRank():
        print(f"[NIRB] Number of snapshot computed : {N}")

    print(Deltas_conv)

    return S, Xi_train, Deltas_conv