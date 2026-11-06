from pathlib import Path

import pytest

import feelpp.core as fppc
import feelpp.core.discr as discr
import feelpp.core.interpolation as interpolation


@pytest.mark.parametrize("family,dim", [("Pch", 1), ("Pch", 2), ("Pchv", 2),
                                      ("Pch", 3), ("Pchv", 3)])
def test_interpolation_order_registrations(family, dim):
    for range_kind in ("range", "range_ptr"):
        classes = [getattr(interpolation,
                           f"OperatorInterpolation_{range_kind}_{family}_{dim}D_P{order}")
                   for order in (1, 2, "Dynamic")]
        assert len(set(classes)) == 3


@pytest.fixture(scope="module")
def interpolation_meshes(init_feelpp, tmp_path_factory):
    fppc.Environment.changeRepository(directory="pyfeelpp-tests/interpolation")
    comm = fppc.Environment.worldComm().to_comm()
    directory = Path(comm.bcast(str(tmp_path_factory.mktemp("interpolation")), root=0))
    meshes = []
    for name, h in (("coarse", 0.4), ("fine", 0.25)):
        geometry = directory / f"{name}.geo"
        if fppc.Environment.isMasterRank():
            geometry.write_text('SetFactory("OpenCASCADE");\n'
                                f'Mesh.MeshSizeMin = {h};\n'
                                f'Mesh.MeshSizeMax = {h};\n'
                                'Rectangle(1) = {0, 0, 0, 1, 2};\n'
                                'Physical Surface("Omega") = {1};\n'
                                'Physical Curve("Boundary") = {1, 2, 3, 4};\n')
        comm.barrier()
        meshes.append(fppc.load(fppc.mesh(dim=2), str(geometry), h))
    assert fppc.nelements(fppc.elements(meshes[0])) != fppc.nelements(fppc.elements(meshes[1]))
    return meshes


@pytest.mark.parametrize("family", ["Pch", "Pchv"])
@pytest.mark.parametrize("dynamic,order", [(False, 1), (False, 2),
                                          (True, 1), (True, 2), (True, 3)])
def test_interpolation_linear_field(interpolation_meshes, family, dynamic, order):
    def make_space(mesh):
        if dynamic:
            return fppc.functionSpace(mesh=mesh, space=family, order=order)
        constructor = getattr(discr, f"{family}_2D_P{order}_G1")
        return constructor(mesh=mesh, worldsComm=fppc.Environment.worldsComm(1))

    domain, image = map(make_space, interpolation_meshes)
    assert domain.order() == image.order() == order
    vectorial = family == "Pchv"
    expression = fppc.expr("{x+y,x-y}:x:y" if vectorial else "x+y:x:y",
                           row=2 if vectorial else 1)
    source = domain.element()
    source.on(range=fppc.elements(domain.mesh()), expr=expression)
    expected = image.element()
    image_range = fppc.elements(image.mesh())
    expected.on(range=image_range, expr=expression)

    operator = interpolation.interpolator(domain=domain, image=image, range=image_range)
    result = operator.interpolate(source)
    assert fppc.normL2(range=image_range, expr=result - expected) < 1e-10
