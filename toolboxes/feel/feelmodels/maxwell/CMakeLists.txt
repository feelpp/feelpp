foreach(_DIM 2 3)
  message(STATUS "[feelpp] MAXWELL ${_DIM}D P1G1")#" N=1..${MAX_ORDER}")
  genLibMaxwell(
    DIM ${_DIM}
    P_ORDER 1 #${_T_ORDER}
    GEO_ORDER 1
    )
  set( FEELPP_MODELS_LIBRARIES_MAXWELL ${FEELPP_MODELS_LIBRARIES_MAXWELL} ${MAXWELL_LIB_NAME} )
endforeach()

add_library( feelpp_toolbox_maxwell_lib SHARED libmaxwell.cpp )
add_library(Feelpp::feelpp_toolbox_maxwell_lib ALIAS feelpp_toolbox_maxwell_lib )  # to match exported target
set_target_properties(feelpp_toolbox_maxwell_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_maxwell_lib")
target_link_libraries( feelpp_toolbox_maxwell_lib ${FEELPP_MODELS_LIBRARIES_MAXWELL} )
set_property( TARGET feelpp_toolbox_maxwell_lib PROPERTY MACOSX_RPATH ON)

foreach(_DIM 2 3)
  message(STATUS "[feelpp] MAXWELL ${_DIM}D P1G1")#" N=1..${MAX_ORDER}")
  genLibMaxwellH1(
    DIM ${_DIM}
    P_ORDER 1
    GEO_ORDER 1
    )
  set( FEELPP_MODELS_LIBRARIES_MAXWELLH1 ${FEELPP_MODELS_LIBRARIES_MAXWELLH1} ${MAXWELLH1_LIB_NAME} )
endforeach()

add_library( feelpp_toolbox_maxwell_h1_lib SHARED libmaxwell_h1.cpp )
add_library(Feelpp::feelpp_toolbox_maxwell_h1_lib ALIAS feelpp_toolbox_maxwell_h1_lib )  # to match exported target
set_target_properties(feelpp_toolbox_maxwell_h1_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_maxwell_h1_lib")
target_link_libraries( feelpp_toolbox_maxwell_h1_lib ${FEELPP_MODELS_LIBRARIES_MAXWELLH1} )
set_property( TARGET feelpp_toolbox_maxwell_h1_lib PROPERTY MACOSX_RPATH ON)

install(TARGETS feelpp_toolbox_maxwell_lib feelpp_toolbox_maxwell_h1_lib DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL( FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/maxwell COMPONENT Devel)
