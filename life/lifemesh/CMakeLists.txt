add_library( lifemesh
  marker.cpp bareitems.cpp 
  entities.cpp simplex.cpp simplexproduct.cpp 
  regiontree.cpp kdtree.cpp
  functors.cpp 
  meshbase.cpp
    )

target_link_libraries (lifemesh lifecore)
