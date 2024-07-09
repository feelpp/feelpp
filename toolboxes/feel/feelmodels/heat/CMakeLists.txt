
feelpp_toolboxes_add_library( heatbase SRCS heatboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )


# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D P<N>G1 N=${_T_ORDER}")
    genLibHeat(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})
  # check that start order is greater than stop order for RANGE below otherwise cmake may have an undefined behavior
  if(${MAX_ORDER} GREATER_EQUAL 2 )
    message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNGN N=2..${MAX_ORDER}")
    foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
      message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D P<N>G<N> N=${_T_ORDER}")
      genLibHeat(
        DIM     ${_DIM}
        T_ORDER ${_T_ORDER}
        GEO_ORDER ${_T_ORDER}
        )
      set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME} )
    endforeach()
  endif()
endforeach()

feelpp_toolboxes_add_library(heat_lib SRCS libheat.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_HEAT})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/heat COMPONENT Devel)
