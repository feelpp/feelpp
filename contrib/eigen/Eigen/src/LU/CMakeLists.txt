FILE(GLOB Eigen_LU_SRCS "*.h")

INSTALL(FILES 
  ${Eigen_LU_SRCS}
  DESTINATION include/feel/Eigen/src/LU COMPONENT Devel
  )

ADD_SUBDIRECTORY(arch)
