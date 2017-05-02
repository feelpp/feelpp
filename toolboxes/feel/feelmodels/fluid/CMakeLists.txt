
# Taylor Hood geo order 1
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FLUID ${_DIM}D PNPN-1G1 N=2..${MAX_ORDER}")
  foreach( _U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _P_ORDER '${_U_ORDER}-1')
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      P_CONTINUITY Continuous
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# Taylor Hood iso-parametric with pressure
foreach( _DIM 2 3)
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FLUID ${_DIM}D PN+1PNGN N=2..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _U_ORDER '${_P_ORDER}+1')
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER ${_P_ORDER}
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      P_CONTINUITY Continuous
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# PnPn geo 1 : same order velocity/pressure
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 1 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FLUID ${_DIM}D PNPNG1 N=1..${MAX_ORDER}")
  foreach( _U_ORDER RANGE 1 ${MAX_ORDER} )
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_U_ORDER}
      P_CONTINUITY Continuous
      ADD_CMAKE_INSTALL
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_fluidmechanics SHARED libfluidmechanics.cpp)
target_link_libraries(feelpp_model_fluidmechanics ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/toolboxes/feel/feelmodels/fluid COMPONENT Devel)