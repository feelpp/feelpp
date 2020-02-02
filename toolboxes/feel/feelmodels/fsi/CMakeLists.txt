


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

add_library(feelpp_toolbox_fsi_lib SHARED libfsi.cpp)
add_library(Feelpp::feelpp_toolbox_fsi_lib ALIAS feelpp_toolbox_fsi_lib )  # to match exported target
set_target_properties(feelpp_toolbox_fsi_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_fsi_lib")
target_link_libraries(feelpp_toolbox_fsi_lib ${FEELPP_MODELS_LIBRARIES_FSI})
set_property(TARGET feelpp_toolbox_fsi_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_fsi_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_fsi_lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/fsi COMPONENT Devel)
