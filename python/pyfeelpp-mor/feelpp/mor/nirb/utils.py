# -*- coding: utf-8 -*-
## utilis function for NIRB
## Thomas Saigre, Ali Elarif
## 09/2022

from errno import EINPROGRESS
import os
import feelpp
import feelpp.mor as mor
import feelpp.toolboxes.core as core
from petsc4py import PETSc
from slepc4py import SLEPc
import numpy as np


############################################################################################################
#                                                                                                          #
# Feel++ utils functions                                                                                   #
#                                                                                                          #
############################################################################################################


def init_feelpp_environment(toolboxType, config_file, argv=['feelpp-mor-nirb']):
    """Initialize Feel++ environment

    Args:
        toolboxType (str): Toolbox used
        config_file (str): path to cfg file
        argv (list, optional): list of arguments to give to Feel++ environment. Defaults to ['feelpp-mor-nirb'].

    Returns:
        feelpp._core.Environment: Environment
    """
    config = feelpp.globalRepository(f"nirb/{toolboxType}")
    opts = core.toolboxes_options(toolboxType).add(mor.makeToolboxMorOptions())
    e = feelpp.Environment(argv, config=config, opts=opts)
    e.setConfigFile(config_file)
    return e



def loadParameterSpace(model_path):
    """Load parameter space from given model path

    Args:
        model_path (str): path to the CRB model file (JSON)

    Returns:
        ParameterSpace: object ParameterSpace
    """

    crb_model_properties = mor.CRBModelProperties("", feelpp.Environment.worldCommPtr(), "")
    crb_model_properties.setup(model_path)
    Dmu = feelpp.mor._mor.ParameterSpace.New(crb_model_properties.parameters(), feelpp.Environment.worldCommPtr())
    return Dmu


def assembleToolbox(tb, mu):
    """Assemble the toolbox tb for the parameter mu

    Args:
        tb (Toolbox): Toolbox object
        mu (parameterSpaceElement): parameter
    """

    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    for i in range(0,mu.size()):
        tb.addParameterInModelProperties(mu.parameterName(i), mu(i))

    tb.updateParameterValues()


def samplingEqui(model_path,Ns,type='equidistribute'):
    """ Get an equidistribute sampling of parameter 

    Args:
        model_path (str): model path 
        Ns (int): number of snapshot 
        type (str, optional): type of the equidistribution. Defaults to 'equidistribute'.
    """
    model = feelpp.readJson(model_path)
    crb = model['CRBParameters']

    Dmu = loadParameterSpace(model_path)
    mu = Dmu.element()

    key = crb.keys()
    mus = []
    dic = {}
    for i in range(Ns): 
        for k in key:
            dic[k] = crb[k]['min'] + i*(crb[k]['max'] - crb[k]['min'])/float(Ns)
        mu.setParameters(dic)
        mus.append(mu)
    
    return mus 


############################################################################################################
#                                                                                                          #
# PETSc handling functions                                                                                 #
#                                                                                                          #
############################################################################################################


def LoadPetscArrayBin(filename):
    """Load a PETSc array from filename, written in binary format

    Args:
        filename (str): path to file to load

    Returns:
        petsc4py.Mat: loaded array
    """
    outputfile = os.path.join(filename)
    viewer = PETSc.Viewer().createBinary(outputfile, 'r')
    # PetscAray = PETSc.Mat().createDense()
    PetscAray = PETSc.Mat().load(viewer)
    return PetscAray

def SavePetscArrayBin(filename, PetscAray):
    """Save a PETSc array to filename, in binary format

    Args:
        filename (str): path to file to save
        PetscAray (petsc4py.Mat): array to save
    """
    outputfile = os.path.join(filename)

    viewer = PETSc.Viewer().createBinary(outputfile, 'w')
    
    # print(help(viewer.pushFormat))
    # print(help(viewer.Format))

    # viewer.pushFormat(viewer.Format.NATIVE)
    # viewer.PetscViewerPushFormat(PETSC_VIEWER_NATIVE)
    viewer(PetscAray)


def TruncatedEigenV(matrix, epsilon = None, nbModes = None):
    """
    Computes a truncated eigen value decomposition of a symetric definite
    matrix in petsc.mat format. Get only eigen value lambda_i > epsilon^2 lambda_max
    the basis vectors returned are orthonormalized 

    Parameters
    ----------
    matrix (petsc.Mat) : the input matrix
    epsilon (float)    : the truncation tolerence, determining the number of keps eigenvalues
    nbModes (int)      : the number of keps eigenvalues

    Returns
    -------
    eigenvalues np.array : kept eigenvalues, of size (nev)
    eigenvectors (list) : kept eigenvectors associated to kept eigenvalues
    """

    if epsilon != None and nbModes != None:# pragma: no cover
        raise("cannot specify both epsilon and nbModes")

    # Get eigenpairs of the matrix 
    E = SLEPc.EPS() # SVD for singular value decomposition or EPS for Eigen Problem Solver  
    E.create()  # create the solver

    E.setOperators(matrix)
    E.setFromOptions()
    E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
    E.setDimensions(matrix.size[1]) # set the number of eigen val to compute
    # E.setTolerances(epsilon) # set the tolerance used for the convergence 

    E.solve()
    nbmaxEv = E.getConverged() # number of eigenpairs 
 
    eigenValues = []
    eigenVectors = E.getInvariantSubspace() # Get orthonormal basis associated to eigenvalues 

    for i in range(nbmaxEv):
        eigenValues.append(float(E.getEigenvalue(i).real))
    
    E.destroy() # destroy the solver object 

    eigenValues = np.array(eigenValues)

    idx = eigenValues.argsort()[::-1]

    eigenValues = eigenValues[idx]
    eigenVectors = [eigenVectors[i] for i in idx]

    # if nbModes == None:
    #     if epsilon == None:
    #         nbModes  = matrix.size[0]
    #     else:
    #         nbModes = 0
    #         bound = (epsilon ** 2) * eigenValues[0]
    #         for e in eigenValues:
    #             if e > bound:
    #                 nbModes += 1
    #         id_max2 = 0
    #         bound = (1 - epsilon ** 2) * np.sum(eigenValues)
    #         temp = 0
    #         for e in eigenValues:
    #             temp += e
    #             if temp < bound:
    #                 id_max2 += 1  # pragma: no cover

    #         nbModes = max(nbModes, id_max2)

    # if nbModes > matrix.size[0]:
    #     print("nbModes taken to max possible value of "+str(matrix.shape[0])+" instead of provided value "+str(nbModes))
    #     nbModes = matrix.size[0]

    # index = np.where(eigenValues<0)
    # if len(eigenValues[index])>0:
    #     if index[0][0]<nbModes:
    #         #print(nbModes, index[0][0])
    #         print("removing numerical noise from eigenvalues, nbModes is set to "+str(index[0][0])+" instead of "+str(nbModes))
    #         nbModes = index[0][0]
    
    return eigenValues, eigenVectors
    # return eigenValues[0:nbModes], eigenVectors[0:nbModes]

############################################################################################################
#                                                                                                          #
# IO handling functions                                                                                    #
#                                                                                                          #
############################################################################################################
def WriteVecAppend(filename, array):
    """ Write an array or list in filename with append mode 
            the vector value will be writen horizontally 
    """
    with open(filename, 'a+') as file:
        file.write(' '.join(str(i) for i in list(array))+"\n")