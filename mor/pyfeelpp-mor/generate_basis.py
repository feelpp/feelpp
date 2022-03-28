import feelpp.mor.reducedbasis.reducedbasis as mor_rb
import sys, os
from feelpp.toolboxes.heat import *
from feelpp.toolboxes.core import *
from feelpp.mor import *
import feelpp

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


import argparse

HOME = os.environ['HOME']

parser = argparse.ArgumentParser()
parser.add_argument('--config-file', help="path to cfg file", type=str)
parser.add_argument('--odir', help="path to output directory", type=str)
parser.add_argument('--dim', help="dimension of the case", type=int)
parser.add_argument('--time-dependant', help="time dependand case", type=bool, default=False)
parser.add_argument('--greedy', help="compute using greedy algorithm", type=bool, default=True)
parser.add_argument('--train-size', help="size of the (random) training set", type=int, default=40)

if __name__ == '__main__':
    args = parser.parse_args()
    dim = args.dim
    config_file = args.config_file
    time_dependant = args.time_dependant
    odir = args.odir
    compute_greedy = args.greedy
    size = args.train_size
else:
    dim = None
    config_file = None
    time_dependant = None
    odir = None
    compute_greedy = None
    size = None


if config_file is not None:
    split = config_file.split('/')
    case = '/'.join(split[:-1])
    if odir is None:
        dir = HOME + '/feel/reducedbasis/'+split[-2]
    else:
        dir = odir
    casefile = split[-1]


def generate_basis():

    if dim not in [2,3]:
        raise ValueError("dim must be 2 or 3")
    if config_file is None:
        raise ValueError("invalid path to config-file")


    if rank == 0:
        print("==============================================")
        print("Generation of the reduced basis")
        print("           Config-file :", f"{case}/{casefile}")
        print("Data will be stored in :", dir)
        print("===============================================")

    config = feelpp.globalRepository(f'{dir}')
    sys.argv = ['generate-basis']
    o = toolboxes_options("heat")
    o.add(makeToolboxMorOptions())

    e = feelpp.Environment(sys.argv, opts=o, config=config)

    feelpp.Environment.setConfigFile(f'{case}/{casefile}')
    feelpp.Environment.changeRepository(directory=f'{dir}/{case}')

    # Set the toolboxes
    heatBox = heat(dim=dim, order=1)
    heatBox.init()



    # model = toolboxmor_2d() if DIM==2 else toolboxmor_3d()
    model = toolboxmor(dim=dim, time_dependent=time_dependant)
    model.setFunctionSpaces( Vh=heatBox.spaceTemperature() )


    def assembleDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleRhs()

    def assembleMDEIM(mu):
        for i in range(0,mu.size()):
            heatBox.addParameterInModelProperties(mu.parameterName(i), mu(i))
        heatBox.updateParameterValues()
        return heatBox.assembleMatrix()

    model.setAssembleDEIM(fct=assembleDEIM)


    model.setAssembleDEIM(fct=assembleDEIM)
    model.setAssembleMDEIM(fct=assembleMDEIM)

    model.initModel()

    heatBoxDEIM = heat(dim=dim, order=1)
    meshDEIM = model.getDEIMReducedMesh()
    heatBoxDEIM.setMesh(meshDEIM)
    heatBoxDEIM.init()

    heatBoxMDEIM = heat(dim=dim, order=1)
    meshMDEIM = model.getMDEIMReducedMesh()
    heatBoxMDEIM.setMesh(meshMDEIM)
    heatBoxMDEIM.init()

    def assembleOnlineDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxDEIM.updateParameterValues()
        return heatBoxDEIM.assembleRhs()

    def assembleOnlineMDEIM(mu):
        for i in range(0,mu.size()):
            heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
        heatBoxMDEIM.updateParameterValues()
        return heatBoxMDEIM.assembleMatrix()

    model.setOnlineAssembleDEIM(assembleOnlineDEIM)
    model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

    model.postInitModel()
    model.setInitialized(True)

    Dmu = model.parameterSpace()


    def listOfParams(n):
        mus = []
        for _ in range(n):
            mus.append(Dmu.element(True, True))
        return mus

    mubar = Dmu.element(True, False)
    # this is temporary : we will be able to read those values from the json
    # mubar.setParameters({"E":40, "T_amb":298, "T_bl":310, "epsilon":1, "h_amb":10, "h_bl":65})
    mubar.setParameters({"Bi":0.01, "k_1":0.1, "k_2":0.1, "k_3":0.1, "k_4":0.1})
    # mubar.setParameters({"phi":250000, "hconv":300, "Tref":450, "P":2e8, "k_4":0.1})
    if rank == 0:
        mubar.view()

    affineDecomposition = model.getAffineDecomposition()
    Aq = affineDecomposition[0]
    Fq = affineDecomposition[1]

    rb = mor_rb.reducedbasis(mor_rb.convertToPetscMat(Aq[0]),
                             mor_rb.convertToPetscVec(Fq[0][0]),
                             model, mubar)
    rb.setVerbose(False)
    if rank == 0:
        print("Size of the big problem :", rb.NN)


    mus = listOfParams(size)
    if compute_greedy:
        mu_0 = Dmu.element()
        rb.greedy(mu_0, mus)
    else:
        rb.computeOfflineReducedBasis(mus)

    rb.computeOfflineErrorRhs()
    rb.computeOfflineError()

    if rank == 0:
        rb.saveReducedBasis(dir, force=True)

    return rb



if __name__ == '__main__':
    generate_basis()