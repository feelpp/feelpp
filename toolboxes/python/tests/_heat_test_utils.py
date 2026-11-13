import feelpp.core as fppc
from pathlib import Path
from feelpp.toolboxes.core import simulate
from feelpp.toolboxes.heat import heat


def run_heat_case(casefile, dim, order, order_geometry=1, repository=None):
    fppc.Environment.setConfigFile(casefile)
    if repository is not None:
        fppc.Environment.changeRepository(directory=repository)
    toolbox = heat(dim=dim, order=order, orderGeometry=order_geometry)
    toolbox.setModelProperties(str(Path(casefile).with_suffix(".json")))
    if not toolbox.isStationary():
        toolbox.setTimeFinal(10 * toolbox.timeStep())
    simulate(toolbox)
    assert toolbox.mesh().order() == order_geometry
    assert toolbox.spaceTemperature().order() == order
    try:
        import pandas as pd

        pd.DataFrame([toolbox.postProcessMeasures().values()])
    except ImportError:
        print("cannot import pandas, no problem it was just a test")
    assert toolbox.checkResults()
