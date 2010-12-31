FILE(GLOB Eigen_LU_arch_SRCS "*.h")

INSTALL(FILES
  ${Eigen_LU_arch_SRCS}
  DESTINATION include/feel/Eigen/src/LU/arch COMPONENT Devel
  )
