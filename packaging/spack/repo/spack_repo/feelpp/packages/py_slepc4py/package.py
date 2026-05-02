from spack.package import *  # noqa: F403
from spack_repo.builtin.packages.py_slepc4py.package import PySlepc4py as BuiltinPySlepc4py


class PySlepc4py(BuiltinPySlepc4py):
    """Feel++ slepc4py overlay fixes for shared Spack environments."""

    # Same class of fix as py-petsc4py's 3.21-3.23 LDSHARED patch: Python's
    # linker command can leak Spack compiler-wrapper paths into the final link.
    patch("ldshared_3211.patch", when="@3.21:3.23")
