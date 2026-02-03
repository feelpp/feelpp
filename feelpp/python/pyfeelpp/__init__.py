import importlib
import sys
import warnings

warnings.warn(
    "The 'pyfeelpp' package is deprecated. Use 'feelpp.core' instead.",
    DeprecationWarning,
    stacklevel=2,
)

core = importlib.import_module("feelpp.core")
sys.modules.setdefault("pyfeelpp.core", core)

try:
    toolboxes = importlib.import_module("feelpp.toolboxes")
except Exception:
    toolboxes = None
else:
    sys.modules.setdefault("pyfeelpp.toolboxes", toolboxes)

try:
    mor = importlib.import_module("feelpp.mor")
except Exception:
    mor = None
else:
    sys.modules.setdefault("pyfeelpp.mor", mor)

from feelpp.core import *

__all__ = getattr(core, "__all__", [])
