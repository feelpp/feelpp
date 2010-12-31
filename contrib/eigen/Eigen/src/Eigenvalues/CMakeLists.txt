FILE(GLOB Eigen_EIGENVALUES_SRCS "*.h")

INSTALL(FILES
  ${Eigen_EIGENVALUES_SRCS}
  DESTINATION include/feel/Eigen/src/Eigenvalues COMPONENT Devel
  )
