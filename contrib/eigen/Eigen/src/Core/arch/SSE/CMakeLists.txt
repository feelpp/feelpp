FILE(GLOB Eigen_Core_arch_SSE_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_arch_SSE_SRCS}
  DESTINATION include/Eigen/src/Core/arch/SSE COMPONENT Devel
)
