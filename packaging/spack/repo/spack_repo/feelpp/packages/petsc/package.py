from spack.package import *  # noqa: F403
from spack_repo.builtin.packages.petsc.package import Petsc as BuiltinPetsc


class Petsc(BuiltinPetsc):
    """Feel++ PETSc overlay fixes for shared Spack environments."""

    def install(self, spec, prefix):
        self.revert_kokkos_nvcc_wrapper()

        # PETSc's lighter install target omits the HPDDM runtime helper library
        # (`libhpddm_petsc`), which PCHPDDM loads dynamically at runtime.
        target = "install" if "+examples" in spec or "+hpddm" in spec else "install-lib"
        make(target, parallel=False)

        if self.run_tests:
            make(f'check PETSC_ARCH="" PETSC_DIR={prefix}', parallel=False)
