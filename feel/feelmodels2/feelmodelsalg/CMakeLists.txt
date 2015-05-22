
add_library( feelmodelsalg SHARED
             modelalgebraic.cpp )
target_link_libraries (feelmodelsalg feelmodelscore ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )
