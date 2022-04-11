
#feelpp_toolboxes_add_library( fsibase SRCS fsiboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

# Taylor Hood geo order 1 , solid : P1G1
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D PNPN-1G1 N=2..${MAX_ORDER}")
  foreach( _FLUID_U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _FLUID_P_ORDER "${_FLUID_U_ORDER}-1")
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER  1
      SOLID_DISP_ORDER 1
      SOLID_GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${FSI_LIB_NAME}  )
  endforeach()
  genLibFSI(
    DIM ${_DIM}
    FLUID_U_ORDER 2
    FLUID_P_ORDER 1
    FLUID_P_CONTINUITY Continuous
    FLUID_GEO_ORDER  1
    SOLID_DISP_ORDER 2
    SOLID_GEO_ORDER 1
    )
  set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${FSI_LIB_NAME}  )
endforeach()

# Taylor Hood iso-parametric with pressure , solid same order pressusre
foreach( _DIM 2 3)
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D iso-parametric PN+1PNGN N=2..${MAX_ORDER}")
  foreach( _FLUID_P_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _FLUID_U_ORDER "${_FLUID_P_ORDER}+1")
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_P_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER ${_FLUID_P_ORDER}
      SOLID_DISP_ORDER ${_FLUID_P_ORDER}
      SOLID_GEO_ORDER ${_FLUID_P_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${FSI_LIB_NAME}  )
  endforeach()
endforeach()

# add Taylor Hood P2P1G2 P2G2
if ( FEELPP_INSTANTIATION_ORDER_MAX GREATER 1 AND FEELPP_MESH_MAX_ORDER GREATER 1 )
  foreach( _DIM 2 3)
    message(STATUS "[feelpp] FSI ${_DIM}D P2P1G2 P2G2")
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER 2
      FLUID_P_ORDER 1
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER 2
      SOLID_DISP_ORDER 2
      SOLID_GEO_ORDER 2
      )
  endforeach()
endif()

# same order velocty,pressure,displacement (geo1)
foreach( _DIM 2 3)
  feelpp_max( MAX_ORDER 1 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] FSI ${_DIM}D same order PNPNPNG1 N=1..${MAX_ORDER}")
  foreach( _FLUID_U_ORDER RANGE 1 ${MAX_ORDER} )
    math(EXPR _FLUID_P_ORDER "${_FLUID_U_ORDER}")
    genLibFSI(
      DIM ${_DIM}
      FLUID_U_ORDER ${_FLUID_U_ORDER}
      FLUID_P_ORDER ${_FLUID_U_ORDER}
      FLUID_P_CONTINUITY Continuous
      FLUID_GEO_ORDER  1
      SOLID_DISP_ORDER ${_FLUID_U_ORDER}
      SOLID_GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_FSI  ${FEELPP_MODELS_LIBRARIES_FSI} ${FSI_LIB_NAME}  )
  endforeach()
endforeach()

feelpp_toolboxes_add_library(fsi_lib SRCS libfsi.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_FSI})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/fsi COMPONENT Devel)
