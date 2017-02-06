
# geo order 1 (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_max( MAX_ORDER 1 ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 1 ${MAX_ORDER} )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER 1
      #PERIODICITY
      POINT_SET EquiSpaced 
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER 1
      #PERIODICITY
      POINT_SET Fekete 
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# isoparametric (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 2 ${MAX_ORDER} )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER ${_POLY_ORDER}
      POINT_SET EquiSpaced 
      #PERIODICITY
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME} )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER ${_POLY_ORDER}
      POINT_SET Fekete 
      #PERIODICITY
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME} )
  endforeach()
endforeach()

# vectorial geo order 1 (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_max( MAX_ORDER 1 ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNVecG1 N=1..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 1 ${MAX_ORDER} )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      POLY_SET Vectorial
      GEO_ORDER 1
      #PERIODICITY
      POINT_SET EquiSpaced 
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      POLY_SET Vectorial
      GEO_ORDER 1
      #PERIODICITY
      POINT_SET Fekete 
      ADD_CMAKE_INSTALL
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_advection SHARED libadvection.cpp )
target_link_libraries(feelpp_model_advection ${FEELPP_MODELS_LIBRARIES_ADVECTION})
