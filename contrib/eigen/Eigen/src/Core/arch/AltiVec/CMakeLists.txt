FILE(GLOB Eigen_Core_arch_AltiVec_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_arch_AltiVec_SRCS}
  DESTINATION include/Eigen/src/Core/arch/AltiVec COMPONENT Devel
)
