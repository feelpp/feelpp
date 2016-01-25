
add_library( feelpp_modelalg SHARED
  modelalgebraicfactory.cpp )
target_link_libraries (feelpp_modelalg ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} feelpp_modelcore )
INSTALL(TARGETS feelpp_modelalg DESTINATION lib/ COMPONENT LibsFeelppModelsCommon)
