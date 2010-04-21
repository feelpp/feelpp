

add_subdirectory (lifecore)
add_subdirectory (lifealg)
add_subdirectory (lifemesh)
#add_subdirectory (lifepoly)
add_subdirectory (lifediscr)
add_subdirectory (lifefilters)
add_subdirectory (lifematerial)
add_subdirectory (lifevf )
add_subdirectory (lifesystem)

add_library(  life SHARED options.cpp )
set_target_properties( life PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})

target_link_libraries (life lifematerial lifevf lifefilters lifediscr lifealg lifemesh lifecore )
set_target_properties(life PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})
