import re

from ._core import Environment
from ._discr import *

_SPACE_ALIASES = {
    "RTh": "Dh",
}

_SPACE_NAME_RE = re.compile(
    r"^(?P<family>[A-Za-z][A-Za-z0-9]*)_(?P<dim>[1-3])D_P(?P<order>Dynamic|[0-9]+)_G(?P<geo>[0-9]+)$"
)


def _discover_space_registry(namespace):
    registry = {}
    legacy = {}
    for name, obj in namespace.items():
        if not isinstance(obj, type):
            continue
        match = _SPACE_NAME_RE.match(name)
        if not match:
            continue
        family = match.group("family")
        dim = int(match.group("dim"))
        geo = int(match.group("geo"))
        order_token = match.group("order")
        entry = registry.setdefault((family, dim, geo), {"dynamic": None, "static": {}})
        if order_token == "Dynamic":
            entry["dynamic"] = obj
            legacy[f"{family}({dim},dynamic,{geo})"] = obj
        else:
            order = int(order_token)
            entry["static"][order] = obj
            legacy[f"{family}({dim},{order},{geo})"] = obj
    return registry, legacy


space_registry, spaces = _discover_space_registry(globals())
space_aliases = dict(_SPACE_ALIASES)


def available_function_spaces():
    return space_registry


def _normalize_order(order):
    try:
        normalized = int(order)
    except (TypeError, ValueError):
        raise RuntimeError(f"FunctionSpace order {order} is invalid, expected integer")
    if normalized < 0:
        raise RuntimeError(f"FunctionSpace order {normalized} is invalid, expected >= 0")
    return normalized


def function_space(mesh, space="Pch", order=1, worldscomm=None):
    """create a function space"""
    if worldscomm is None:
        worldscomm = Environment.worldsComm(1)

    order = _normalize_order(order)
    family = _SPACE_ALIASES.get(space, space)
    dim = mesh.dimension()
    geo = mesh.order()
    entry = space_registry.get((family, dim, geo))

    if entry is not None:
        static_ctor = entry["static"].get(order)
        # Keep dynamic spaces as default, except Pdh(order=0) where legacy APIs
        # (eg pid()) require the static P0 type.
        if family == "Pdh" and order == 0 and static_ctor is not None:
            return static_ctor(mesh=mesh, worldsComm=worldscomm, runtimeOrder=order)
        if entry["dynamic"] is not None:
            return entry["dynamic"](mesh=mesh, worldsComm=worldscomm, runtimeOrder=order)
        if static_ctor is not None:
            return static_ctor(mesh=mesh, worldsComm=worldscomm, runtimeOrder=order)
        available_orders = [str(k) for k in sorted(entry["static"].keys())]
        if entry["dynamic"] is not None:
            available_orders.insert(0, "dynamic")
        raise RuntimeError(
            f"FunctionSpace {family}({dim},{order},{geo}) not available, "
            f"supported orders: {', '.join(available_orders)}"
        )

    key = f"{family}({dim},{order},{geo})"
    if key in spaces:
        return spaces[key](mesh=mesh, worldsComm=worldscomm, runtimeOrder=order)
    dkey = f"{family}({dim},dynamic,{geo})"
    if dkey in spaces:
        return spaces[dkey](mesh=mesh, worldsComm=worldscomm, runtimeOrder=order)

    available_families = sorted({k[0] for k in space_registry.keys()})
    raise RuntimeError(
        f"FunctionSpace {family}({dim},{order},{geo}) not available. "
        f"Known families: {', '.join(available_families)}"
    )


# Backward-compatible aliases (legacy camelCase API)
def availableFunctionSpaces():
    return available_function_spaces()


def functionSpace(mesh, space="Pch", order=1, worldscomm=None):
    return function_space(mesh=mesh, space=space, order=order, worldscomm=worldscomm)
