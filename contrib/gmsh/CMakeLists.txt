message(STATUS "[feelpp] adding missing Gmsh headers...")

# we include this directory : add some missing headers from Gmsh
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/contrib/gmsh )

install(FILES BasisFactory.h OpenFile.h DESTINATION include/feel/gmsh )
