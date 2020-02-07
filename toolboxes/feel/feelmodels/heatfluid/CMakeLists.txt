
foreach(_DIM 2 3)
  feelpp_max( MAX_ORDER 2 ${FEELPP_INSTANTIATION_ORDER_MAX})
  message(STATUS "[feelpp] HEATFLUID ${_DIM}D PNPN-1P1G1 N=2..${MAX_ORDER}")
  message(STATUS "[feelpp] HEATFLUID ${_DIM}D PNPN-1PNG1 N=2..${MAX_ORDER}")
  # geo order 1, Taylor-Hood, T_ORDER=U_ORDER and T_ORDER=1
  foreach( _U_ORDER RANGE 2 ${MAX_ORDER} )
    math(EXPR _P_ORDER "${_U_ORDER}-1")
    genLibHeatFluid(
      DIM     ${_DIM}
      T_ORDER 1
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEATFLUID ${FEELPP_MODELS_LIBRARIES_HEATFLUID} ${HEATFLUID_LIB_NAME}  )
    genLibHeatFluid(
      DIM     ${_DIM}
      T_ORDER ${_U_ORDER}
      U_ORDER ${_U_ORDER}
      P_ORDER ${_P_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEATFLUID ${FEELPP_MODELS_LIBRARIES_HEATFLUID} ${HEATFLUID_LIB_NAME}  )
  endforeach()
  message(STATUS "[feelpp] HEATFLUID ${_DIM}D PNPNPNG1 N=1..1")
  # geo order 1, T_ORDER=U_ORDER=P_ORDER
  foreach( _U_ORDER RANGE 1 1 )
    genLibHeatFluid(
      DIM     ${_DIM}
      T_ORDER ${_U_ORDER}
      U_ORDER ${_U_ORDER}
      P_ORDER ${_U_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_HEATFLUID ${FEELPP_MODELS_LIBRARIES_HEATFLUID} ${HEATFLUID_LIB_NAME}  )
  endforeach()
endforeach()


add_library(feelpp_toolbox_heatfluid_lib SHARED libheatfluid.cpp )
add_library(Feelpp::feelpp_toolbox_heatfluid_lib ALIAS feelpp_toolbox_heatfluid_lib )  # to match exported target
set_target_properties(feelpp_toolbox_heatfluid_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_heatfluid_lib")
target_link_libraries(feelpp_toolbox_heatfluid_lib ${FEELPP_MODELS_LIBRARIES_HEATFLUID})
set_property(TARGET feelpp_toolbox_heatfluid_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_heatfluid_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_heatfluid_lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/heatfluid COMPONENT Devel)
