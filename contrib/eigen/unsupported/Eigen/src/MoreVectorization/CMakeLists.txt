FILE(GLOB Eigen_MoreVectorization_SRCS "*.h")

INSTALL(FILES
  ${Eigen_MoreVectorization_SRCS}
  DESTINATION include/feel/unsupported/Eigen/src/MoreVectorization COMPONENT Devel
  )
