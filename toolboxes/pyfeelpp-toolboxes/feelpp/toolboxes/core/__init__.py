from ._modelcore import *
from ._modelmesh import *
import feelpp


def simulate(toolbox, export=True, buildModelAlgebraicFactory=True,data=None):
    """simulate a toolbox

    simulate execute the toolbox in steady or transient case and export the results

    Parameters:
        toolbox -- a toolbox which has been configured

    Example::
        feelpp.Environment.setConfigFile('toolboxs/laplace/l-shape/l-shape.cfg')
        toolbox = toolboxs.toolboxs(dim=2)
        simulate(toolbox)        
    """
    toolbox.init(buildModelAlgebraicFactory)
    #toolbox.printAndSaveInfo()
    if toolbox.isStationary():
        toolbox.solve()
        if export:
            toolbox.exportResults()
    else:
        if not toolbox.doRestart():
            toolbox.exportResults(toolbox.timeInitial())
        toolbox.startTimeStep()
        while not toolbox.timeStepBase().isFinished():
            if feelpp.Environment.isMasterRank():
                print("============================================================\n")
                print("time simulation: {}s/{}s with step: {}".format(toolbox.time(),toolbox.timeFinal(),toolbox.timeStep()))
                print("============================================================\n")
            toolbox.solve()
            if export:
                toolbox.exportResults()
            toolbox.updateTimeStep()
    return not toolbox.checkResults()
