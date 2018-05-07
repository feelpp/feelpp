math(EXPR _T_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
feelpp_max( MAX_ORDER 1 ${_T_ORDER})
message(STATUS "[feelpp] MAXWELL 3D PNG1 N=1..${MAX_ORDER}")
foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
  genLibMaxwell(
    DIM 3
    P_ORDER ${_T_ORDER}
    GEO_ORDER 1
    )
  set( FEELPP_MODELS_LIBRARIES_MAXWELL ${FEELPP_MODELS_LIBRARIES_MAXWELL} ${MAXWELL_LIB_NAME} )
endforeach()

message(STATUS "maxwell lib: '${FEELPP_MODELS_LIBRARIES_MAXWELL}'" )

add_library( feelpp_toolbox_maxwell_lib SHARED libmaxwell.cpp )
target_link_libraries( feelpp_toolbox_maxwell_lib ${FEELPP_MODELS_LIBRARIES_MAXWELL} )
set_property( TARGET feelpp_toolbox_maxwell_lib PROPERTY MACOSX_RPATH ON)
install(TARGETS feelpp_toolbox_maxwell_lib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL( FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/maxwell COMPONENT Devel)
