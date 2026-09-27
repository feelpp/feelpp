from spack_repo.builtin.packages.mpich.package import Mpich as BuiltinMpich


class Mpich(BuiltinMpich):
    """Use LUMI's documented CH3 configure flag for MPICH 3.4.3."""

    def configure_args(self):
        args = super().configure_args()
        if self.spec.satisfies("@3.4.3 device=ch3 netmod=tcp"):
            return ["--with-device=ch3" if arg == "--with-device=ch3:nemesis:tcp" else arg for arg in args]
        return args
