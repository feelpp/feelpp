FILE(GLOB Eigen_QR_SRCS "*.h")

INSTALL(FILES
  ${Eigen_QR_SRCS}
  DESTINATION include/Eigen/src/QR COMPONENT Devel
  )
