import feelpp.core as fppc
from feelpp.toolboxes.cfpdes import cfpdes
from feelpp.toolboxes.core import simulate
from _case_paths import toolbox_case


def run_cfpdes_case(casefile, dim):
    fppc.Environment.setConfigFile(toolbox_case(casefile))
    toolbox = cfpdes(dim=dim)
    if not toolbox.isStationary():
        toolbox.setTimeFinal(toolbox.timeStep() * 10)
    simulate(toolbox)
    assert toolbox.checkResults()
