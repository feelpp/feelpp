# This file generate a bash script for the feel++ interpreter.
# It configure feel++ libraries for cling, but not a builtin cling executable!
# NB: This file is intended to be included in a custom target.

# Get a list of subdirectory.
macro(subdirlist result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            list(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
endmacro()

#subdirlist( FEELPP_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/feel )
feelpp_expand_target_libraries( FEELPP_ALL ${FEELPP_LIBRARIES} )

# Create the script feel++.

# Script head.
file( WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"#!/usr/bin/env bash
printf \"\\\n\\e[1m\\e[1;31m\
 _____ _____ _____ __      _     _
|   __|   __|   __|  |   _| |_ _| |
|   __|   __|   __|  |__|_   _|_   _|
|__|  |_____|_____|_____| |_|   |_|
\\e[0m \\e[1m\\e[1;36m
Feel++ interpreter (cling based)
For more infos, see www.feelpp.org
Type '.help' for help, '.q' to exit
\\e[0m
\"
"
)

# Script cling binary program.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"${CLING_BIN} \\"
)

# Set environment variables.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-DBOOST_FILESYSTEM_VERSION=3 \\
-DBOOST_NO_CXX11_SCOPED_ENUMS \\
-DBOOST_NO_SCOPED_ENUMS \\
-DBOOST_OPTIONAL_USE_OLD_DEFINITION_OF_NONE=1 \\
-DBOOST_PARAMETER_MAX_ARITY=24 \\
-DBOOST_PP_VARIADICS=0 \\
-DBOOST_RESULT_OF_USE_TR1 \\
-DBOOST_TEST_DYN_LINK \\
-DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR \\
-DFEELPP_HAS_DLFCN_H \\
-DFEELPP_HAS_DLOPEN \\
-DFEELPP_HAS_GMSH=1 \\
-DFEELPP_HAS_GMSH_ADAPT_H \\
-DFEELPP_HAS_GSL=1 \\
-DFEELPP_HAS_HPDDM \\
-DFEELPP_HAS_JSONLAB \\
-DFEELPP_HAS_METIS \\
-DFEELPP_HAS_MPI=1 \\
-DFEELPP_HAS_MPI_H=1 \\
-DFEELPP_HAS_PARALUTION \\
-DFEELPP_HAS_PETSC \\
-DFEELPP_HAS_PETSC_H \\
-DGFLAGS_IS_A_DLL=0 \\
-DGOOGLE_GLOG_DLL_DECL=\"\" \\
-DGOOGLE_GLOG_DLL_DECL_FOR_UNITTESTS=\"\" \\
-DHAVE_LIBDL \\
-DIN_GINAC \\
"
)

# Include precompiled headers.
if( FEELPP_ENABLE_PCH )
    file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-include-pch ${CMAKE_BINARY_DIR}/feel/cotire/feelpp_CXX_prefix.hxx.pch \\
"
#-fsyntax-only \\
#-femit-all-decls \\
    )
endif()

# Include path.
set( FEELPP_INTERPRETER_INCLUDE_DIRS
    ${FEELPP_INCLUDE_DIR}
    ${FEELPP_DEPS_INCLUDE_DIR}
    # contrib source directories.
    ${google-glog_SOURCE_DIR}/src
    ${gflags_SOURCE_DIR}
    ${GiNaC_SOURCE_DIR}
    ${GiNaC_SOURCE_DIR}/ginac
    ${CMAKE_SOURCE_DIR}/contrib/nlopt/api
    ${Eigen3_SOURCE_DIR}
    ${Eigen3_SOURCE_DIR}/unsupported
    ${METIS_SOURCE_DIR}
    # contrib binary directories.
    ${google-glog_BINARY_DIR}
    ${gflags_BINARY_DIR}/include
    ${GiNaC_BINARY_DIR}
    ${GiNaC_BINARY_DIR}/ginac
    ${Eigen3_BINARY_DIR}
    ${Eigen3_BINARY_DIR}/unsupported
    ${METIS_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/contrib/nlopt/api
)
list( REMOVE_DUPLICATES FEELPP_INTERPRETER_INCLUDE_DIRS )

# Library dirs path or libraries path (eg. -L/path/to/lib.so).
set( FEELPP_INTERPRETER_LIBRARY_DIRS
       ${CMAKE_INSTALL_PREFIX}/lib
       ${CMAKE_BINARY_DIR}/feel
    #   ${FEELPP_DEPS_LINK_DIR}
    # Contrib binary directories.
    #   ${google-glog_BINARY_DIR}
    #   ${gflags_BINARY_DIR}
    #   ${GiNaC_BINARY_DIR}
    #   ${GiNaC_BINARY_DIR}/ginac
    #   ${Eigen3_BINARY_DIR}
    #   ${METIS_BINARY_DIR}/libmetis
    #   ${CMAKE_BINARY_DIR}/contrib/nlopt/
    # Variable generated from feelpp_find_library
    #   ${FEELPP_ALL_LIBRARY_DIRS}
)
list( REMOVE_DUPLICATES FEELPP_INTERPRETER_LIBRARY_DIRS )


# Libraries to link (eg. -ltoto).
set( FEELPP_INTERPRETER_LINK_LIBRARIES
    #${FEELPP_LINK_LIBRARIES}
    ${FEELPP_LIBRARY} # -lfeepp
)
list( REMOVE_DUPLICATES FEELPP_INTERPRETER_LINK_LIBRARIES )


# Shared libraries to link (eg. /path/to/toto.so)
set( FEELPP_INTERPRETER_LIBRARIES
    "" #empty
    #${FEELPP_LIBRARIES}
)
list( REMOVE_DUPLICATES FEELPP_INTERPRETER_LIBRARIES )


foreach( incdir ${FEELPP_INTERPRETER_INCLUDE_DIRS} )
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-I${incdir} \\
"
)
endforeach()

foreach( libdir ${FEELPP_INTERPRETER_LIBRARY_DIRS} )
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-L${libdir} \\
"
)
endforeach()

foreach( lib ${FEELPP_INTERPRETER_LINK_LIBRARIES} )
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-l${lib} \\
"
)
endforeach()

foreach( lib ${FEELPP_INTERPRETER_LIBRARIES} )
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
${lib} \\
"
)
endforeach()

# Script cling options.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"--std=c++${FEELPP_STD_CPP} \\
-ftemplate-depth=1024 \\
--nologo \\
$1"
)

# Set permissions as executable script.
file( INSTALL ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    DESTINATION ${CLING_INSTALL_PREFIX}
    )
