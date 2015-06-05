
# Taylor Hood geo order 1 , solid : P1G1
foreach( _DIM 2 3)
  foreach( _FLUID_U_ORDER 2 3 4 5 )
    math(EXPR _FLUID_P_ORDER '${_FLUID_U_ORDER}-1')
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER  1
      SOLID_DISP_ORDER 1
      SOLID_GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# Taylor Hood iso-parametric with pressure , solid same order pressusre
foreach( _DIM 2 3)
  foreach( _FLUID_U_ORDER 3 4 5 )
    math(EXPR _FLUID_P_ORDER '${_FLUID_U_ORDER}-1')
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER ${_FLUID_P_ORDER}
      SOLID_DISP_ORDER ${_FLUID_P_ORDER}
      SOLID_GEO_ORDER ${_FLUID_P_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_fsi SHARED libfsi.cpp)
target_link_libraries(feelpp_model_fsi ${FEELPP_MODELS_LIBRARIES_FSI})
