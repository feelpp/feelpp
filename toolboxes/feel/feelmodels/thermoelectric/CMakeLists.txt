
# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
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
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
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

add_library(feelpp_model_thermoelectric SHARED libthermoelectric.cpp )
target_link_libraries(feelpp_model_thermoelectric ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC})
set_property(TARGET feelpp_model_thermoelectric PROPERTY MACOSX_RPATH ON)
install(TARGETS feelpp_model_thermoelectric DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/thermoelectric COMPONENT Devel)
