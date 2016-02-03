


# Taylor Hood geo order 1 , solid : P1G1
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D PN+1PNG1 N=2..${MAX_ORDER}")
  foreach( _FLUID_U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _FLUID_P_ORDER '${_FLUID_U_ORDER}-1')
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER  1
      SOLID_DISP_ORDER 1
      SOLID_GEO_ORDER 1
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# Taylor Hood iso-parametric with pressure , solid same order pressusre
foreach( _DIM 2 3)
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D iso-parametric PN+1PNGN N=2..${MAX_ORDER}")
  foreach( _FLUID_P_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _FLUID_U_ORDER '${_FLUID_P_ORDER}+1')
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER ${_FLUID_P_ORDER}
      SOLID_DISP_ORDER ${_FLUID_P_ORDER}
      SOLID_GEO_ORDER ${_FLUID_P_ORDER}
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# same order velocty,pressure,displacement (geo1)
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 1 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D same order PNPNPNG1 N=1..${MAX_ORDER}")
  foreach( _FLUID_U_ORDER RANGE 1 ${MAX_ORDER} )
    math(EXPR _FLUID_P_ORDER '${_FLUID_U_ORDER}')
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_U_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER  1
      SOLID_DISP_ORDER ${_FLUID_U_ORDER}
      SOLID_GEO_ORDER 1
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_fsi SHARED libfsi.cpp)
target_link_libraries(feelpp_model_fsi ${FEELPP_MODELS_LIBRARIES_FSI})
