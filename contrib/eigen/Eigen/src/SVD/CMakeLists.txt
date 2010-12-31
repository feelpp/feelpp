FILE(GLOB Eigen_SVD_SRCS "*.h")

INSTALL(FILES
  ${Eigen_SVD_SRCS}
  DESTINATION include/feel/Eigen/src/SVD COMPONENT Devel
  )
