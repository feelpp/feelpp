FILE(GLOB Eigen_Jacobi_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Jacobi_SRCS}
  DESTINATION include/feel/Eigen/src/Jacobi COMPONENT Devel
  )
