FILE(GLOB Eigen_plugins_SRCS "*.h")

INSTALL(FILES
  ${Eigen_plugins_SRCS}
  DESTINATION include/Eigen/src/plugins COMPONENT Devel
  )
