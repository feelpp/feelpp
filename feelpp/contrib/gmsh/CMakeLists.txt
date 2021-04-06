
# add some missing headers from Gmsh
if( FEELPP_HAS_GMSH_LIBRARY AND ( NOT FEELPP_HAS_GMSH_API OR FEELPP_GMSH_API_MISSING_HEADERS ) )
  message(STATUS "[feelpp] adding missing Gmsh headers...")

  if ( FEELPP_GMSH_API_MISSING_HEADERS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}.${GMSH_PATCH_VERSION} )
    message(STATUS "[feelpp] adding missing Gmsh headers from contrib")
    set(FEELPP_GMSH_MISSING_HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}.${GMSH_PATCH_VERSION} PARENT_SCOPE)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}.${GMSH_PATCH_VERSION}/gmsh.h DESTINATION include/feelpp/gmsh COMPONENT Devel )
  else()

    find_program(WGET NAMES wget )
    if ( WGET )
      message(STATUS "[feelpp] adding missing Gmsh headers from gitlab downloads")
      #INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )
      #foreach( _F Numeric/BasisFactory.h CommonOpenFile.h Mesh/Field.h )
      add_custom_target(feelpp_gmsh_get_files ALL)
      #add_dependencies(feelpp_contrib feelpp_gmsh_get_files)
      set(gmsh_gitlab_url "http://gitlab.onelab.info/gmsh/gmsh/raw/gmsh_${GMSH_MAJOR_VERSION}_${GMSH_MINOR_VERSION}_${GMSH_PATCH_VERSION}")
      set(FEELPP_GMSH_MISSING_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)
      if ( NOT FEELPP_HAS_GMSH_API )
        include(CheckIncludeFileCXX)
        CHECK_INCLUDE_FILE_CXX("gmsh/Field.h" FEELPP_HAS_GMSH_FIELD_H)
        CHECK_INCLUDE_FILE_CXX("gmsh/BasisFactory.h" FEELPP_HAS_GMSH_BASISFACTORY_H)
        CHECK_INCLUDE_FILE_CXX("gmsh/OpenFile.h" FEELPP_HAS_GMSH_OPENFILE_H)
        message(STATUS "[feelpp] gmsh headers xz${FEELPP_HAS_GMSH_FIELD_H} ${FEELPP_HAS_GMSH_BASISFACTORY_H} ${FEELPP_HAS_GMSH_OPENFILE_H}")
        if ( NOT FEELPP_HAS_GMSH_FIELD_H OR NOT FEELPP_HAS_GMSH_BASISFACTORY_H OR NOT FEELPP_HAS_GMSH_OPENFILE_H )
          add_custom_command(
            TARGET feelpp_gmsh_get_files
            COMMAND ${WGET} -c -q --timeout=5 ${gmsh_gitlab_url}/Mesh/Field.h -O Field.h 
            COMMAND ${WGET} -c -q --timeout=5 ${gmsh_gitlab_url}/Numeric/BasisFactory.h -O BasisFactory.h
            COMMAND ${WGET} -c -q --timeout=5 ${gmsh_gitlab_url}/Common/OpenFile.h -O OpenFile.h )
          install(FILES ${CMAKE_CURRENT_BINARY_DIR}/BasisFactory.h ${CMAKE_CURRENT_BINARY_DIR}/Field.h ${CMAKE_CURRENT_BINARY_DIR}/OpenFile.h DESTINATION include/feelpp/gmsh COMPONENT Devel )
        endif()
      else()
        #file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        add_custom_command(
          TARGET feelpp_gmsh_get_files
          COMMAND ${WGET} -c  --timeout=10 ${gmsh_gitlab_url}/api/gmsh.h -O gmsh.h || true )
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gmsh.h DESTINATION include/feelpp/gmsh COMPONENT Devel )
      endif()
    endif()
  endif()
endif()

