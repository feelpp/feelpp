FILE(GLOB Eigen_Householder_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Householder_SRCS}
  DESTINATION include/Eigen/src/Householder COMPONENT Devel
  )
