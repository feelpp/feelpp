import feelpp.core as fppc
from feelpp.toolboxes.core import simulate
from feelpp.toolboxes.heat import heat


def run_heat_case(casefile, dim, order):
    fppc.Environment.setConfigFile(casefile)
    toolbox = heat(dim=dim, order=order)
    if not toolbox.isStationary():
        toolbox.setTimeFinal(10 * toolbox.timeStep())
    simulate(toolbox)
    try:
        import pandas as pd

        pd.DataFrame([toolbox.postProcessMeasures().values()])
    except ImportError:
        print("cannot import pandas, no problem it was just a test")
    assert toolbox.checkResults()
