
add_library( feelpp_modelalg SHARED
             modelalgebraic.cpp )
target_link_libraries (feelpp_modelalg ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} feelpp_modelcore )
