
feelpp_toolboxes_add_library( heatbase SRCS heatboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )


# Keep only P1/G1 at compile time; all other orders are runtime values.
foreach(_DIM 2 3)
  foreach(_ORDER 1 Dynamic)
    message(STATUS "[feelpp] HEAT ${_DIM}D P${_ORDER}G${_ORDER}")
    genLibHeat(
      DIM ${_DIM}
      T_ORDER ${_ORDER}
      GEO_ORDER ${_ORDER}
      )
    list(APPEND FEELPP_MODELS_LIBRARIES_HEAT ${HEAT_LIB_NAME})
  endforeach()
endforeach()

feelpp_toolboxes_add_library(heat_lib SRCS libheat.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_HEAT})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/heat COMPONENT Devel)
