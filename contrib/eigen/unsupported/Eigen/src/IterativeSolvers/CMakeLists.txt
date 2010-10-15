FILE(GLOB Eigen_IterativeSolvers_SRCS "*.h")

INSTALL(FILES
  ${Eigen_IterativeSolvers_SRCS}
  DESTINATION include/unsupported/Eigen/src/IterativeSolvers COMPONENT Devel
  )
