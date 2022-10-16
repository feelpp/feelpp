from ._modelcore import *
from ._modelmesh import *
import feelpp


def simulate(toolbox, export=True, buildModelAlgebraicFactory=True,data=None):
    """simulate a toolbox

    simulate execute the toolbox in steady or transient case and export the results

    Parameters:
        toolbox -- a toolbox which has been configured

    Returns:
        a tuple (status,results) where status is a boolean(True if success, False otherwise) and results is a dictionary

    Example::
        feelpp.Environment.setConfigFile('toolboxs/laplace/l-shape/l-shape.cfg')
        toolbox = toolboxs.toolboxs(dim=2)
        [success,measures]=simulate(toolbox)        
        import pandas as pd
        df = pd.DataFrame(measures)
        print(df.head())
    """
    toolbox.init(buildModelAlgebraicFactory)
    #toolbox.printAndSaveInfo()
    meas=[]
    if toolbox.isStationary():
        toolbox.solve()
        if export:
            toolbox.exportResults()
        meas.append(toolbox.postProcessMeasures().values())
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
            meas.append(toolbox.postProcessMeasures().values())
            if export:
                toolbox.exportResults()
            toolbox.updateTimeStep()
    return [toolbox.checkResults(),meas]
