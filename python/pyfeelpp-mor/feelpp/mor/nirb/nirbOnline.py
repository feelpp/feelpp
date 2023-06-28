from time import time
import feelpp
from feelpp.mor.nirb import nirbOnline
from feelpp.mor.nirb.utils import WriteVecAppend, init_feelpp_environment
import pandas as pd
import os
from feelpp.mor.nirb.nirb_perf import *
import argparse
from feelpp.timing import tic,toc


def run_online(nirb_on: nirbOnline, Nmu=10, Xi=None, Nb=-1, export=False, computeError=False, rectification=True):
    """Run NIRB online step

    Parameters
    ----------
    nirb_on : nirbOnline
        nirbOnline object, with the basis initialized
    Nmu : int, optional
        number of parameters if none is given, by default 10
    Xi : list of parameterSpaceElement, optional
        list of parameters used to compute, by default None. If None, a sampling of size Nmu is generated
    Nb : int, optional
        size of the basis, by default -1. If -1, the whole basis is used
    export : bool, optional
        export solutions for visualization, by default False
    computeError : bool, optional
        compute error, by default False
    rectification : bool, optional
        use rectification post process, by default True
    """
    
    if Xi is None:
        s = nirb_on.Dmu.sampling()
        s.sampling(Nmu, "log-random")
        Xi = s.getVector()
    
    wo = ['without', 'with'][rectification]

    if export:
        dirname = "nirbOnlineSolutions"
        nirb_on.initExporter(dirname, toolbox="fine")

    for i, mu in enumerate(Xi):
        tic()
        uHhN = nirb_on.getOnlineSolution(mu, Nb, doRectification=rectification)
        toc(f"NIRB : getOnlineSolution {wo} rectification")
        if feelpp.Environment.isMasterRank():
            print(f"[NIRB Online] Solution computed for mu = {mu}, with a basis of size {Nb}")

        if export:
            nirb_on.exportField(uHhN, f"n{i:02}uHhN{mu}")

        if computeError:
            uh = nirb_on.getToolboxSolution(nirb_on.tbFine, mu)
            diff = uHhN - uh
            error = nirb_on.normMat(diff)
            print(f"[NIRB Online] L2 norm between nirb online and toolbox sol = {error}")
            if export:
                nirb_on.exportField(diff, f"n{i:02}diff_{mu}")
                nirb_on.exportField(uh, f"n{i:02}uh{mu}")

    if export:
        print(f"[NIRB Online] Exporting nirb solutions in {dirname}")
        nirb_on.saveExporter()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='NIRB Online')
    parser.add_argument('--config-file', type=str, help='path to cfg file')
    parser.add_argument("--Nsnap", help="Size of basis to load [default=-1]. If -1, the whole basis is loaded", type=int, default=-1)
    parser.add_argument("--N", help="Size of basis use for online computations [default=-1]. If -1, the whole basis is used", type=int, default=-1)

    parser.add_argument("--rectification", help="With or without rectification [default=1]", type=int, default=1)
    parser.add_argument("--outdir", help="output directory", type=str, default=None)
    parser.add_argument("--export", help="Export nirb sol for vizualisation [default=0]", type=int, default=0)
    parser.add_argument("--compute-error", help="Compute errors of nirb solutions [default=0]", type=int, default=0)

    parser.add_argument("--Xi", help="Path to file containing the parameters", type=str, default=None)
    parser.add_argument("--Nparam", help="Number of parameters to test, if Xi is not given [default=10]", type=int, default=10)

    args = parser.parse_args()
    config_file = args.config_file
    outdir = args.outdir

    cfg = feelpp.readCfg(config_file)
    toolboxType = cfg['nirb']['toolboxType']
    e = init_feelpp_environment(toolboxType, config_file)

    nirb_file = feelpp.Environment.expand(cfg['nirb']['filename'])
    config_nirb = feelpp.readJson(nirb_file)['nirb']


    convergence = args.convergence != 0
    export = args.export != 0
    computeError = args.compute_error != 0
    doRectification = args.rectification != 0

    nbSnap = args.Nsnap
    Nb = args.N

    doGreedy = config_nirb['doGreedy']
    doRectification = config_nirb['doRectification']
    rectPath = ["noRect", "Rect"][doRectification]
    greedyPath = ["noGreedy", "Greedy"][doGreedy]

    if outdir is None:
        RESPATH = f"results/{rectPath}/{greedyPath}"
    else:
        RESPATH = outdir

    start = time()
    nirb_on = nirbOnline(**config_nirb)
    nirb_on.initModel()
    err = nirb_on.loadData(path=RESPATH, nbSnap=nbSnap)
    finish = time()
    assert err == 0, "Error while loading data"

    if feelpp.Environment.isMasterRank():
        print(f"[NIRB Online] Data loaded in {finish-start} s")

    s = nirb_on.Dmu.sampling()
    Xi = None
    if args.Xi is not None:
        N_train = s.readFromFile('./sampling.sample')
        Xi = s.getVector()

    run_online(nirb_on, Nmu=args.Nparam, Xi=Xi, Nb=Nb, export=export, computeError=computeError, rectification=doRectification)


    feelpp.Environment.saveTimers(display=True)
