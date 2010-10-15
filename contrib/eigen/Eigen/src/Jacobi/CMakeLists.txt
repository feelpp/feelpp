FILE(GLOB Eigen_Jacobi_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Jacobi_SRCS}
  DESTINATION include/Eigen/src/Jacobi COMPONENT Devel
  )
