FILE(GLOB Eigen_AutoDiff_SRCS "*.h")

INSTALL(FILES
  ${Eigen_AutoDiff_SRCS}
  DESTINATION include/unsupported/Eigen/src/AutoDiff COMPONENT Devel
  )
