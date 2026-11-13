
feelpp_toolboxes_add_library( electricbase SRCS electricboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# Keep only P1/G1 at compile time; all other orders are runtime values.
foreach(_DIM 2 3)
  foreach(_ORDER 1 Dynamic)
    message(STATUS "[feelpp] ELECTRIC ${_DIM}D P${_ORDER}G${_ORDER}")
    genLibElectric(
      DIM ${_DIM}
      P_ORDER ${_ORDER}
      GEO_ORDER ${_ORDER}
      )
    list(APPEND FEELPP_MODELS_LIBRARIES_ELECTRIC ${ELECTRIC_LIB_NAME})
  endforeach()
endforeach()

feelpp_toolboxes_add_library(electric_lib SRCS libelectric.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_ELECTRIC})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/electric COMPONENT Devel)
