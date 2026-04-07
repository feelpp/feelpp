import feelpp.core as fppc
import feelpp.core.quality as fppcq
from feelpp.toolboxes.cfpdes import cfpdes
from feelpp.toolboxes.core import simulate
from _case_paths import toolbox_case


def test_cfpdes_remesh():
    fppc.Environment.changeRepository(directory="pyfeelpptoolboxes-tests/cfpdes/laplace/l-shape")
    fppc.Environment.setConfigFile(toolbox_case("cfpdes/laplace/l-shape/l-shape.cfg"))
    toolbox = cfpdes(dim=2)
    simulate(toolbox, export=False)

    exporter = fppc.exporter(mesh=toolbox.mesh(), name="l-shape", geo="change")
    exporter.step(0.0).setMesh(toolbox.mesh())
    toolbox.exportSolutionToStep(exporter.step(0.0))

    hclose = 0.05
    hfar = 1

    function_space = fppc.functionSpace(mesh=toolbox.mesh())
    metric = fppc.gradedls(function_space, fppc.boundaryfaces(function_space.mesh()), hclose, hfar)
    exporter.step(0.0).add("metric", metric)
    exporter.step(0.0).add("quality", fppcq.etaQ(toolbox.mesh()))
    exporter.save()

    new_mesh, _ = fppc.remesh(
        mesh=toolbox.mesh(),
        metric="gradedls({},{})".format(hclose, hfar),
        required_elts=[],
        required_facets=[],
        parent=None,
    )

    remeshed_toolbox = cfpdes(dim=2)
    remeshed_toolbox.setMesh(new_mesh)
    simulate(remeshed_toolbox, export=False)
    exporter.step(1.0).setMesh(new_mesh)
    remeshed_toolbox.exportSolutionToStep(exporter.step(1.0))

    remeshed_function_space = fppc.functionSpace(mesh=remeshed_toolbox.mesh())
    remeshed_metric = fppc.gradedls(
        remeshed_function_space,
        fppc.boundaryfaces(remeshed_function_space.mesh()),
        hclose,
        hfar,
    )
    quality = fppcq.etaQ(remeshed_toolbox.mesh())
    exporter.step(1.0).add("metric", remeshed_metric)
    exporter.step(1.0).add("quality", quality)
    exporter.save()
    assert remeshed_toolbox.checkResults()
