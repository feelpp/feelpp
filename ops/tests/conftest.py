from __future__ import annotations

from pathlib import Path
import os
import sys


OPS_SRC = Path(__file__).resolve().parents[1] / "src"
REPO_ROOT = Path(__file__).resolve().parents[2]

if str(OPS_SRC) not in sys.path:
    sys.path.insert(0, str(OPS_SRC))

os.environ.setdefault("FEELPP_REPO_ROOT", str(REPO_ROOT))

