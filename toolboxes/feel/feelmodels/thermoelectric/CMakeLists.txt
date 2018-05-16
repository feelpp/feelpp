
# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibThermoelectric(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      P_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})
  message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
    genLibThermoelectric(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      P_ORDER ${_T_ORDER}
      GEO_ORDER ${_T_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME} )
  endforeach()
endforeach()

if( ${FEELPP_MESH_MAX_ORDER} GREATER 1 )
  foreach(_DIM 2 3)
    math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
    feelpp_max( MAX_ORDER 1 ${_T_ORDER})
    message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNG2 N=2..${MAX_ORDER}")
    foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
      genLibThermoelectric(
        DIM     ${_DIM}
        T_ORDER ${_T_ORDER}
        P_ORDER ${_T_ORDER}
        GEO_ORDER 2
        )
      set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME} )
    endforeach()
  endforeach()
endif()


add_library(feelpp_toolbox_thermoelectric_lib SHARED libthermoelectric.cpp )
target_link_libraries(feelpp_toolbox_thermoelectric_lib ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC})
set_property(TARGET feelpp_toolbox_thermoelectric_lib PROPERTY MACOSX_RPATH ON)
install(TARGETS feelpp_toolbox_thermoelectric_lib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/thermoelectric COMPONENT Devel)
