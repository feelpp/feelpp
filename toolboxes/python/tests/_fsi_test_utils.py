import feelpp.core as fppc
from feelpp.toolboxes.fsi import fsi
from _case_paths import toolbox_case


def run_fsi_case(casefile, dim, order_u, order_p, order_geo):
    fppc.Environment.setConfigFile(toolbox_case(casefile))
    toolbox = fsi(dim=dim, orderU=order_u, orderP=order_p, orderGeo=order_geo)
    toolbox.init()

    assert toolbox.modelFluid() is not None
    assert toolbox.modelSolid() is not None
    assert toolbox.fieldVelocity() is not None
    assert toolbox.fieldPressure() is not None
    assert toolbox.fieldDisplacement() is not None

    if toolbox.isStationary():
        toolbox.solve()
    else:
        toolbox.setTimeFinal(toolbox.timeInitial() + toolbox.timeStep())
        toolbox.startTimeStep()
        toolbox.solve()
        toolbox.updateTimeStep()

    assert toolbox.checkResults()
