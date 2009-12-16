add_library( lifemesh SHARED
  marker.cpp bareitems.cpp
  entities.cpp simplex.cpp simplexproduct.cpp
  regiontree.cpp kdtree.cpp
  functors.cpp
  meshbase.cpp
    )

target_link_libraries (lifemesh lifecore)
set_target_properties(lifemesh PROPERTIES VERSION ${LIFE_SHARED_VERSION} SOVERSION ${LIFE_SHARED_SOVERSION})
