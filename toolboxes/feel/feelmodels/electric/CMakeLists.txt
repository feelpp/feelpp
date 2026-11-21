
feelpp_toolboxes_add_library( electricbase SRCS electricboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# Static baseline + dynamic runtime-order instantiation
foreach(_DIM 2 3)
  message(STATUS "[feelpp] ELECTRIC ${_DIM}D static P1G1")
  genLibElectric(
    DIM ${_DIM}
    P_ORDER 1
    GEO_ORDER 1
    )
  list(APPEND FEELPP_MODELS_LIBRARIES_ELECTRIC ${ELECTRIC_LIB_NAME})

  message(STATUS "[feelpp] ELECTRIC ${_DIM}D dynamic P/G")
  genLibElectric(
    DIM ${_DIM}
    P_ORDER Dynamic
    GEO_ORDER 1
    )
  list(APPEND FEELPP_MODELS_LIBRARIES_ELECTRIC ${ELECTRIC_LIB_NAME})
endforeach()

feelpp_toolboxes_add_library(electric_lib SRCS libelectric.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_ELECTRIC})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/electric COMPONENT Devel)
