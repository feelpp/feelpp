
add_library( feelpp_modelalg SHARED
  modelalgebraicfactory.cpp )
target_link_libraries (feelpp_modelalg ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} feelpp_modelcore )
set_property(TARGET feelpp_modelalg PROPERTY MACOSX_RPATH ON)
if( FEELPP_ENABLE_PCH_MODELS )
    add_precompiled_header( feelpp_modelalg )
endif()
INSTALL(TARGETS feelpp_modelalg DESTINATION lib/ COMPONENT Libs EXPORT feelpp-toolboxes-targets)
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/modelalg COMPONENT Devel)
