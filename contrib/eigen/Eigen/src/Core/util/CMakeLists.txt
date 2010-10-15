FILE(GLOB Eigen_Core_util_SRCS "*.h")

INSTALL(FILES 
  ${Eigen_Core_util_SRCS}
  DESTINATION include/Eigen/src/Core/util COMPONENT Devel
  )
