FILE(GLOB Eigen_NonLinear_SRCS "*.h")

INSTALL(FILES
  ${Eigen_NonLinear_SRCS}
  DESTINATION include/unsupported/Eigen/src/NonLinear COMPONENT Devel
  )
