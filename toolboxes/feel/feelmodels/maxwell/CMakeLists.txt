foreach(_DIM 2 3)
  # math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  # feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] MAXWELL ${_DIM}D P1G1")#" N=1..${MAX_ORDER}")
  # foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibMaxwell(
      DIM ${_DIM}
      P_ORDER 1 #${_T_ORDER}
      GEO_ORDER 1
      )
    set( FEELPP_MODELS_LIBRARIES_MAXWELL ${FEELPP_MODELS_LIBRARIES_MAXWELL} ${MAXWELL_LIB_NAME} )
  # endforeach()
endforeach()

# if( ${FEELPP_MESH_MAX_ORDER} GREATER 1 )
#   foreach(_DIM 2 3)
#     math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
#     feelpp_max( MAX_ORDER 1 ${_T_ORDER})
#     message(STATUS "[feelpp] MAXWELL ${_DIM}D PNG2 N=2..${MAX_ORDER}")
#     foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
#       genLibMaxwell(
#         DIM     ${_DIM}
#         P_ORDER ${_T_ORDER}
#         GEO_ORDER 2
#         )
#       set(FEELPP_MODELS_LIBRARIES_MAXWELL ${FEELPP_MODELS_LIBRARIES_MAXWELL} ${MAXWELL_LIB_NAME} )
#     endforeach()
#   endforeach()
# endif()

add_library( feelpp_toolbox_maxwell_lib SHARED libmaxwell.cpp )
add_library(Feelpp::feelpp_toolbox_maxwell_lib ALIAS feelpp_toolbox_maxwell_lib )  # to match exported target
set_target_properties(feelpp_toolbox_maxwell_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_maxwell_lib")
target_link_libraries( feelpp_toolbox_maxwell_lib ${FEELPP_MODELS_LIBRARIES_MAXWELL} )
set_property( TARGET feelpp_toolbox_maxwell_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_maxwell_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_maxwell_lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL( FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/maxwell COMPONENT Devel)
