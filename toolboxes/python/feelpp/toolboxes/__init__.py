import importlib

_SUBMODULES = {
    "core": "core",
    "cfpdes": "cfpdes",
    "heat": "heat",
    "electric": "electric",
    "fluid": "fluid",
    "solid": "solid",
    "hdg": "hdg",
    "heatfluid": "heatfluid",
    "thermoelectric": "thermoelectric",
    "fsi": "fsi",
    "advection": "advection",
    "maxwell": "maxwell",
}

__all__ = sorted(_SUBMODULES)


def __getattr__(name):
    if name not in _SUBMODULES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return importlib.import_module(f".{_SUBMODULES[name]}", __name__)


def __dir__():
    return sorted(set(globals()) | set(__all__))
