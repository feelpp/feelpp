import feelpp.core as fppc
from feelpp.toolboxes.solid import solid
from _case_paths import toolbox_case


def run_solid_case(casefile, dim, order_disp):
    fppc.Environment.setConfigFile(toolbox_case(casefile))
    toolbox = solid(dim=dim, orderDisp=order_disp)
    toolbox.init()
    toolbox.solve()
    assert toolbox.mesh() is not None
