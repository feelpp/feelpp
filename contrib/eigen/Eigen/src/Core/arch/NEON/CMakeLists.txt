FILE(GLOB Eigen_Core_arch_NEON_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_arch_NEON_SRCS}
  DESTINATION include/Eigen/src/Core/arch/NEON COMPONENT Devel
)
