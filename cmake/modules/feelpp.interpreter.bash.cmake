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

subdirlist( FEELPP_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/feel )

# Create the script feel++.

# Script head.
file( WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"#!/bin/sh
printf \"\n
 _____ _____ _____ __      _     _
|   __|   __|   __|  |   _| |_ _| |
|   __|   __|   __|  |__|_   _|_   _|
|__|  |_____|_____|_____| |_|   |_|

Feel++ interpreter (cling based)
For more infos, see www.feelpp.org
Type '.help' for help, '.q' to exit

\"
"
)

# Script cling binary program.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"${CLING_BIN} \\"
)

# Script install include directories.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"
-I${GMSH_INCLUDE_DIR} \\
-I${PARALUTION_INCLUDE_DIR} \\
-I${MPI_CXX_INCLUDE_PATH} \\
-I${HPDDM_INCLUDE_DIR} \\
-I${CMAKE_SOURCE_DIR} \\
-I${google-glog_SOURCE_DIR}/src \\
-I${gflags_SOURCE_DIR} \\
-I${GiNaC_SOURCE_DIR} \\
-I${GiNaC_SOURCE_DIR}/ginac \\
-I${CMAKE_SOURCE_DIR}/contrib/nlopt/api \\
-I${Eigen3_SOURCE_DIR} \\
-I${Eigen3_SOURCE_DIR}/unsupported \\
-I${METIS_SOURCE_DIR} \\
-I${PETSC_DIR}/include \\
-I${SLEPC_DIR}/include \\
-I${CMAKE_INSTALL_PREFIX}/include/ \\
"
)

# Script install include directories.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-I${google-glog_BINARY_DIR} \\
-I${gflags_BINARY_DIR}/include \\
-I${GiNaC_BINARY_DIR} \\
-I${GiNaC_BINARY_DIR}/ginac \\
-I${Eigen3_BINARY_DIR} \\
-I${Eigen3_BINARY_DIR}/unsupported \\
-I${METIS_BINARY_DIR} \\
-I${CMAKE_BINARY_DIR}/contrib/nlopt/api \\
"
)

# Script lib directories.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-L${CMAKE_SOURCE_DIR}/feel \\
-L${google-glog_BINARY_DIR} \\
-L${gflags_BINARY_DIR} \\
-L${GiNaC_BINARY_DIR} \\
-L${Eigen3_BINARY_DIR} \\
-L${METIS_BINARY_DIR}/libmetis \\
-L${CMAKE_BINARY_DIR}/contrib/nlopt/ \\
-L${PETSC_DIR}/lib \\
-L${SLEPC_DIR}/lib \\
-L${CMAKE_INSTALL_PREFIX}/lib/ \\
"
)

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
-DGOOGLE_GLOG_DLL_DECL="" \\
-DGOOGLE_GLOG_DLL_DECL_FOR_UNITTESTS="" \\
-DHAVE_LIBDL \\
-DIN_GINAC \\
"
)

## Script source include directories
#foreach( subdir ${FEELPP_INCLUDE_DIRS} )
#    file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
#        "-I${CMAKE_SOURCE_DIR}/${subdir} \\\n"
#        )
#endforeach()

if( FEELPP_ENABLE_PCH )
# Script cling binary program.
    file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"\
-include-pch ${CMAKE_SOURCE_DIR}/feel/feel.hpp.pch \\
"
    )
endif()

# Script cling options.
file( APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
"--std=c++${FEELPP_STD_CPP} \\
-ftemplate-depth=1024 \\
-stdlib=libstdc++ \\
-lfeelpp_metis \\
-lfeelpp_nlopt \\
-lfeelpp_ginac \\
-lfeelpp_glog \\
-lfeelpp \\
--nologo \\
$1"
)

# Set permissions as executable script.
file( INSTALL ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/feel++
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    DESTINATION ${CLING_INSTALL_PREFIX}
    )
