
feelpp_toolboxes_add_library( fluidbase SRCS fluidmechanicsboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# Taylor Hood geo order 1
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FLUID ${_DIM}D PNPN-1G1 N=2..${MAX_ORDER}")
  foreach( _U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _P_ORDER "${_U_ORDER}-1")
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      P_CONTINUITY Continuous
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${FLUIDMECHANICS_LIB_NAME}  )
  endforeach()
endforeach()

# Taylor Hood iso-parametric with pressure
foreach( _DIM 2 3)
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${FEELPP_INSTANTIATION_ORDER_MAX})

  # check that start order is greater than stop order for RANGE below otherwise cmake may have an undefined behavior
  if( ${MAX_ORDER} GREATER_EQUAL 2)
    message(STATUS "[feelpp] FLUID ${_DIM}D PN+1PNGN N=2..${MAX_ORDER}")
    foreach( _P_ORDER RANGE 2 ${MAX_ORDER} )
      math(EXPR _U_ORDER "${_P_ORDER}+1")
      genLibFluidMechanics(
        DIM ${_DIM}
        GEO_ORDER ${_P_ORDER}
        U_ORDER ${_U_ORDER}
        P_ORDER ${_P_ORDER}
        P_CONTINUITY Continuous
        )
      set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${FLUIDMECHANICS_LIB_NAME}  )
    endforeach()
  endif()
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
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${FLUIDMECHANICS_LIB_NAME}  )
  endforeach()
endforeach()

# Taylor Hood P2P1G2
foreach( _DIM 2 3)
  if ( FEELPP_INSTANTIATION_ORDER_MAX GREATER 1 AND FEELPP_MESH_MAX_ORDER GREATER 1 )
    genLibFluidMechanics(
      DIM ${_DIM}
      GEO_ORDER 2
      U_ORDER 2
      P_ORDER 1
      P_CONTINUITY Continuous
      )
    set(FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS  ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS} ${FLUIDMECHANICS_LIB_NAME}  )
  endif()
endforeach()

feelpp_toolboxes_add_library(fluid_lib SRCS libfluidmechanics.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_FLUIDMECHANICS})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/fluid COMPONENT Devel)
