FILE(GLOB Eigen_AutoDiff_SRCS "*.h")

INSTALL(FILES
  ${Eigen_AutoDiff_SRCS}
  DESTINATION include/feel/unsupported/Eigen/src/AutoDiff COMPONENT Devel
  )
