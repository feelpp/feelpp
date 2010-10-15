FILE(GLOB Eigen_BVH_SRCS "*.h")

INSTALL(FILES
  ${Eigen_BVH_SRCS}
  DESTINATION include/unsupported/Eigen/src/BVH COMPONENT Devel
  )
