
# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibHeat(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})
  message(STATUS "[feelpp] HEAT TRANSFER ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
    genLibHeat(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      GEO_ORDER ${_T_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_HEAT  ${FEELPP_MODELS_LIBRARIES_HEAT} ${HEAT_LIB_NAME} )
  endforeach()
endforeach()

add_library(feelpp_toolbox_heat_lib SHARED libheat.cpp )
target_link_libraries(feelpp_toolbox_heat_lib ${FEELPP_MODELS_LIBRARIES_HEAT})
set_property(TARGET feelpp_toolbox_heat_lib PROPERTY MACOSX_RPATH ON)
install(TARGETS feelpp_toolbox_heat_lib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/heat COMPONENT Devel)
