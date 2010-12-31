FILE(GLOB Eigen_Eigen2Support_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Eigen2Support_SRCS}
  DESTINATION include/feel/Eigen/src/Eigen2Support COMPONENT Devel
  )
