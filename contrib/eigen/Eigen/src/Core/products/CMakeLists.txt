FILE(GLOB Eigen_Core_Product_SRCS "*.h")

INSTALL(FILES
  ${Eigen_Core_Product_SRCS}
  DESTINATION include/feel/Eigen/src/Core/products COMPONENT Devel
  )
