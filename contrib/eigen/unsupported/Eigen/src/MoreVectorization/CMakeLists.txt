FILE(GLOB Eigen_MoreVectorization_SRCS "*.h")

INSTALL(FILES
  ${Eigen_MoreVectorization_SRCS}
  DESTINATION include/unsupported/Eigen/src/MoreVectorization COMPONENT Devel
  )
