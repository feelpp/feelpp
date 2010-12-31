FILE(GLOB Eigen_Geometry_arch_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Geometry_arch_SRCS}
  DESTINATION include/feel/Eigen/src/Geometry/arch COMPONENT Devel
  )
