"""Expose Debian's real PETSc/SLEPc bindings in runtime images."""

from pathlib import Path
import sysconfig


multiarch = sysconfig.get_config_var("MULTIARCH")
paths = []
for library in ("petsc", "slepc"):
    root = Path("/usr/lib") / f"{library}dir"
    pattern = f"{library}*/{multiarch}-real/lib/python3/dist-packages/{library}4py/__init__.py"
    candidates = {module.parent.parent.resolve() for module in root.glob(pattern)}
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one installed {library}4py real binding, found {sorted(candidates)}")
    paths.append(str(candidates.pop()))

# Load after the distribution's .pth files so PETSC_DIR/SLEPC_DIR and default
# alternatives keep precedence when development packages provide them.
destination = Path("/usr/lib/python3/dist-packages/zz-feelpp-bindings.pth")
destination.write_text("\n".join(paths) + "\n", encoding="utf-8")
