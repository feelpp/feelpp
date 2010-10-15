FILE(GLOB Eigen_Polynomials_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Polynomials_SRCS}
  DESTINATION include/unsupported/Eigen/src/Polynomials COMPONENT Devel
  )
