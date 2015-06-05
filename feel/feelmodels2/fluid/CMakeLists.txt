
# Taylor Hood geo order 1
foreach( _DIM 2 3)
  foreach( _U_ORDER 2 3 4 5 )
    math(EXPR _P_ORDER '${_U_ORDER}-1')
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      P_CONTINUITY Continuous
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# Taylor Hood iso-parametric with pressure
foreach( _DIM 2 3)
  foreach( _U_ORDER 3 4 5 )
    math(EXPR _P_ORDER '${_U_ORDER}-1')
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER ${_P_ORDER}
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      P_CONTINUITY Continuous
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# PnPn geo 1 : same order velocity/pressure
foreach( _DIM 2 3)
  foreach( _U_ORDER 1 2 )
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_U_ORDER}
      P_CONTINUITY Continuous
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_fluidmechanics SHARED libfluidmechanics.cpp)
target_link_libraries(feelpp_model_fluidmechanics ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS})
