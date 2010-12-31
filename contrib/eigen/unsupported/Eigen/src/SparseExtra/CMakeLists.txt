FILE(GLOB Eigen_SparseExtra_SRCS "*.h")

INSTALL(FILES
  ${Eigen_SparseExtra_SRCS}
  DESTINATION include/feel/unsupported/Eigen/src/SparseExtra COMPONENT Devel
  )
