FILE(GLOB Eigen_SparseExtra_SRCS "*.h")

INSTALL(FILES
  ${Eigen_SparseExtra_SRCS}
  DESTINATION include/unsupported/Eigen/src/SparseExtra COMPONENT Devel
  )
