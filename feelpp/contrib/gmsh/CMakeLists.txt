message(STATUS "[feelpp] adding missing Gmsh headers...")

# we include this directory : add some missing headers from Gmsh
#INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/contrib/gmsh )
find_program(WGET NAMES wget )
if ( WGET )
  #INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )
  #foreach( _F Numeric/BasisFactory.h CommonOpenFile.h Mesh/Field.h )
  add_custom_target(gmsh_get_files ALL)
  set(gmsh_gitlab_url "http://gitlab.onelab.info/gmsh/gmsh/raw/gmsh_${GMSH_MAJOR_VERSION}_${GMSH_MINOR_VERSION}_${GMSH_PATCH_VERSION}")
  add_custom_command(
    TARGET gmsh_get_files
    COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/Field.h -O Field.h 
    COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Numeric/BasisFactory.h -O BasisFactory.h
    COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Common/OpenFile.h -O OpenFile.h )
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/BasisFactory.h ${CMAKE_CURRENT_BINARY_DIR}/Field.h ${CMAKE_CURRENT_BINARY_DIR}/OpenFile.h DESTINATION include/feelpp/gmsh COMPONENT Devel )
endif()

