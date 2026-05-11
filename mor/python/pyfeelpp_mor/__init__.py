import importlib
import warnings

warnings.warn(
    "The 'pyfeelpp_mor' package is deprecated. Use 'feelpp.mor' instead.",
    DeprecationWarning,
    stacklevel=2,
)

_mor = importlib.import_module("feelpp.mor")


def __getattr__(name):
    return getattr(_mor, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_mor)))


__all__ = getattr(_mor, "__all__", [])
