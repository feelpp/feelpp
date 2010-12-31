FILE(GLOB Eigen_StlSupport_SRCS "*.h")

INSTALL(FILES
  ${Eigen_StlSupport_SRCS}
  DESTINATION include/feel/Eigen/src/StlSupport COMPONENT Devel
  )
