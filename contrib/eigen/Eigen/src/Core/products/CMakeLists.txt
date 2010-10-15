FILE(GLOB Eigen_Core_Product_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_Product_SRCS}
  DESTINATION include/Eigen/src/Core/products COMPONENT Devel
  )
