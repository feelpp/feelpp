FILE(GLOB Eigen_Core_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_SRCS}
  DESTINATION include/Eigen/src/Core COMPONENT Devel
  )

ADD_SUBDIRECTORY(products)
ADD_SUBDIRECTORY(util)
ADD_SUBDIRECTORY(arch)
