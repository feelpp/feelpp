FILE(GLOB Eigen_Core_arch_Default_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_arch_Default_SRCS}
  DESTINATION include/Eigen/src/Core/arch/Default COMPONENT Devel
)
