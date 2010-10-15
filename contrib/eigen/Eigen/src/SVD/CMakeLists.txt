FILE(GLOB Eigen_SVD_SRCS "*.h")

INSTALL(FILES
  ${Eigen_SVD_SRCS}
  DESTINATION include/Eigen/src/SVD COMPONENT Devel
  )
