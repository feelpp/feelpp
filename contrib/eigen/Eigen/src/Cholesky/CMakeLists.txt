FILE(GLOB Eigen_Cholesky_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Cholesky_SRCS}
  DESTINATION include/Eigen/src/Cholesky COMPONENT Devel
  )
