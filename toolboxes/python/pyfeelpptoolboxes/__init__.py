import importlib
import warnings

warnings.warn(
    "The 'pyfeelpptoolboxes' package is deprecated. Use 'feelpp.toolboxes' instead.",
    DeprecationWarning,
    stacklevel=2,
)

_toolboxes = importlib.import_module("feelpp.toolboxes")


def __getattr__(name):
    return getattr(_toolboxes, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_toolboxes)))


__all__ = getattr(_toolboxes, "__all__", [])
