FILE(GLOB Eigen_MatrixFunctions_SRCS "*.h")

INSTALL(FILES
  ${Eigen_MatrixFunctions_SRCS}
  DESTINATION include/feel/unsupported/Eigen/src/MatrixFunctions COMPONENT Devel
  )
