FILE(GLOB Eigen_misc_SRCS "*.h")

INSTALL(FILES
  ${Eigen_misc_SRCS}
  DESTINATION include/feel/Eigen/src/misc COMPONENT Devel
  )
