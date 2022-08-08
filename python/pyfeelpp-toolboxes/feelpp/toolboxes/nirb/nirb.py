import feelpp
import feelpp.mor as mor
from feelpp.toolboxes.heat import *
import feelpp.interpolation as fi
import json5 as json



def loadModel(model_path):
    """Load the model from given modle path

    Args:
        model_path (str): path to the model file (JSON)

    Returns:
        json: model loaded
    """
    f = open(model_path, "r")
    model = json.load(f)
    f.close()
    return model


def setToolbox(h, geo_path, model, dim=2, order=2):
    """Set up the toolbox object for the given model and mesh

    Args:
        h (float): size of the mesh
        geo_path (str): path to the goemetries file
        model (str): path to the model file
        dim (int): dimension of the mesh
        order (int): order of the finite elements

    Returns:
        Toolbox: toolbox object
    """

    # load meshes
    mesh_ = feelpp.mesh(dim=2, realdim=2)
    mesh = feelpp.load(mesh_, geo_path, h)

    # set mesh and model properties
    tb = heat(dim=dim, order=order)
    tb.setMesh(mesh)
    tb.setModelProperties(model)

    tb.init()

    return tb

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



def createInterpolator(image_tb,domain_tb):
    """Create an interpolator between two toolboxes
    
    Args:
        source (Toolbox): coarse toolbox
        image (Toolbox): fine toolbox

    Returns:
        OperatorInterpolation: interpolator object
    """
    Vh_image = image_tb.spaceTemperature()
    Vh_domain = domain_tb.spaceTemperature()
    interpolator = fi.interpolator(domain = Vh_domain, image = Vh_image, range = image_tb.rangeMeshElements())
    return interpolator

