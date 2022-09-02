import feelpp
import feelpp.mor as mor
import feelpp.toolboxes.heat as heat
import feelpp.toolboxes.fluid as fluid
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


def setToolbox(h, geo_path, model, dim=2, order=2, type_tb="heat"):
    """Set up the toolbox object for the given model and mesh

    Args:
        h (float): size of the mesh
        geo_path (str): path to the goemetries file
        model (str): path to the model file
        dim (int): dimension of the mesh
        order (int): order of the finite elements
        type_tb (str): name of the toolbox {"heat"|"fluid"}

    Returns:
        Toolbox: toolbox object
    """

    # load meshes
    mesh_ = feelpp.mesh(dim=2, realdim=2)
    mesh = feelpp.load(mesh_, geo_path, h)

    # set mesh and model properties
    if type_tb == "heat":
        tb = heat.heat(dim=dim, order=order)
    elif type_tb == "fluid":
        tb = fluid.fluid(dim=dim)
    else:
        raise ValueError("Unknown toolbox")

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

def getField(toolbox, type_tb):
    """Get field of interest from the toolbox

    Args:
        toolbox (Toolbox): Tolbox object
        type_tb (str): name of the toolbox {"heat"|"fluid"}

    Raises:
        ValueError: Unknow toolbox

    Returns:
        feelpp_.discr.Element_*: field of the solution
    """
    if type_tb == "heat":
        return toolbox.fieldTemperature()
    elif type_tb == "fluid":
        return toolbox.fieldVelocity()
    else:
        raise ValueError("Unknown toolbox")


def createInterpolator(domain_tb, image_tb, type_tb):
    """Create an interpolator between two toolboxes
    
    Args:
        domain_tb (Toolbox): coarse toolbox
        image_tb  (Toolbox): fine toolbox

    Returns:
        OperatorInterpolation: interpolator object
    """
    if type_tb == "heat":
        Vh_image = image_tb.spaceTemperature()
        Vh_domain = domain_tb.spaceTemperature()
    elif type_tb == "fluid":
        Vh_image = image_tb.functionSpaceVelocity()
        Vh_domain = domain_tb.functionSpaceVelocity()
    interpolator = fi.interpolator(domain = Vh_domain, image = Vh_image, range = image_tb.rangeMeshElements())
    return interpolator
