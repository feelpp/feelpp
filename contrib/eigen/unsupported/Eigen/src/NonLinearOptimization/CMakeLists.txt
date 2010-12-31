FILE(GLOB Eigen_NonLinear_SRCS "*.h")

INSTALL(FILES
  ${Eigen_NonLinear_SRCS}
  DESTINATION include/feel/unsupported/Eigen/src/NonLinear COMPONENT Devel
  )
