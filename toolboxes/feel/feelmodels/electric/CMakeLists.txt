
feelpp_toolboxes_add_library( electricbase SRCS electricboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] ELECTRIC ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibElectric(
      DIM     ${_DIM}
      P_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_ELECTRIC  ${FEELPP_MODELS_LIBRARIES_ELECTRIC} ${ELECTRIC_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}-1")
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})

  # check that start order is greater than stop order for RANGE below otherwise cmake may have an undefined behavior
  if(${MAX_ORDER} GREATER_EQUAL 2)
    message(STATUS "[feelpp] ELECTRIC ${_DIM}D PNGN N=2..${MAX_ORDER}")
    foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
      genLibElectric(
        DIM     ${_DIM}
        P_ORDER ${_T_ORDER}
        GEO_ORDER ${_T_ORDER}
        )
      set(FEELPP_MODELS_LIBRARIES_ELECTRIC  ${FEELPP_MODELS_LIBRARIES_ELECTRIC} ${ELECTRIC_LIB_NAME} )
    endforeach()
  endif()
endforeach()

feelpp_toolboxes_add_library(electric_lib SRCS libelectric.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_ELECTRIC})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/electric COMPONENT Devel)
