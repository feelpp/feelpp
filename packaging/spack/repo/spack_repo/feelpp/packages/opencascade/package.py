from spack.package import *  # noqa: F403
from spack_repo.builtin.packages.opencascade.package import Opencascade as BuiltinOpencascade


class Opencascade(BuiltinOpencascade):
    """Feel++ OpenCASCADE overlay fixes for shared Spack environments."""

    def cmake_args(self):
        args = super().cmake_args()

        if self.spec.satisfies("platform=darwin"):
            tcl_libs = self.spec["tcl"].libs
            if tcl_libs:
                args.append(self.define("3RDPARTY_TCL_LIBRARY", tcl_libs[0]))
            tcl_lib_dirs = self.spec["tcl"].libs.directories
            if tcl_lib_dirs:
                args.append(self.define("3RDPARTY_TCL_LIBRARY_DIR", tcl_lib_dirs[0]))

        return args
