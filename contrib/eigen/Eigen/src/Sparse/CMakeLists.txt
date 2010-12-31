FILE(GLOB Eigen_Sparse_SRCS "*.h")

INSTALL(FILES 
  ${Eigen_Sparse_SRCS}
  DESTINATION include/feel/Eigen/src/Sparse COMPONENT Devel
  )
