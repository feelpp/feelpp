
# geo order 1
foreach(_DIM 2 3)
  foreach( _T_ORDER 1 2 3 4 )
    genLibThermoDynamics(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_THERMODYNAMICS  ${FEELPP_MODELS_LIBRARIES_THERMODYNAMICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  foreach( _T_ORDER 2 3 4 )
    genLibThermoDynamics(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER ${_T_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_THERMODYNAMICS  ${FEELPP_MODELS_LIBRARIES_THERMODYNAMICS} ${LIBBASE_NAME} )
  endforeach()
endforeach()

add_library(feelpp_model_thermodynamics SHARED libthermodynamics.cpp )
target_link_libraries(feelpp_model_thermodynamics ${FEELPP_MODELS_LIBRARIES_THERMODYNAMICS})
