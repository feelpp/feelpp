####################################################################
#
# Usage:
#  - create a new folder, let's call it cdash
#  - in that folder, do:
#    ctest -S path/to/feel/cmake/dashboard/testsuite.cmake[,option1=value1[,option2=value2]]
#
# Options:
#  - FEELPP_CXX: compiler, eg.: g++-4.2
#      default: default c++ compiler
#  - FEELPP_SITE: eg, sd-25154, or the name of the contributor, etc.
#      default: hostname
#  - FEELPP_BUILD_STRING: a string which identify the system/compiler. It should be formed like that:
#        <OS_name>-<OS_version>-<arch>-<compiler-version>
#      with:
#        <OS_name> = opensuse, debian, osx, windows, cygwin, freebsd, solaris, etc.
#        <OS_version> = 11.1, XP, vista, leopard, etc.
#        <arch> = i386, x86_64, ia64, powerpc, etc.
#        <compiler-version> = gcc-4.3.2, icc-11.0, MSVC-2008, etc.
#  - FEELPP_CMAKE_DIR: path to cmake executable
#  - FEELPP_MODE: dashboard model, can be Experimental, Nightly, or Continuous
#      default: Nightly
#  - FEELPP_WORK_DIR: directory used to download the source files and make the builds
#      default: folder which contains this script
#  - FEELPP_CMAKE_ARGS: additional arguments passed to cmake
#  - FEELPP_GENERATOR_TYPE: allows to overwrite the generator type
#      default: nmake (windows
#      See http://www.cmake.org/cmake/help/cmake2.6docs.html#section_Generators for a complete
#      list of supported generators.
#  - FEELPP_NO_UPDATE: allows to submit dash boards from local repositories
#      This might be interesting in case you want to submit dashboards
#      including local changes.
#  - CTEST_SOURCE_DIRECTORY: path to feel's feel (use a new and empty folder, not the one you are working on)
#      default: <FEELPP_WORK_DIR>/feel
#  - CTEST_BINARY_DIRECTORY: build directory
#      default: <FEELPP_WORK_DIR>/nightly-<FEELPP_CXX>
#
# Here is an example running several compilers on a linux system:
# #!/bin/bash
# ARCH=`uname -m`
# SITE=`hostname`
# VERSION=opensuse-11.1
# WORK_DIR=/home/prudhomm/scratch/cdash
# # get the last version of the script
# wget http://bitbucket.org/feel/feel/raw/tip/test/testsuite.cmake -o $WORK_DIR/testsuite.cmake
# COMMON="ctest -S $WORK_DIR/testsuite.cmake,FEELPP_WORK_DIR=$WORK_DIR,FEELPP_SITE=$SITE,FEELPP_MODE=$1,FEELPP_BUILD_STRING=$OS_VERSION-$ARCH"
# $COMMON-gcc-3.4.6,FEELPP_CXX=g++-3.4
# $COMMON-gcc-4.0.1,FEELPP_CXX=g++-4.0.1
# $COMMON-icc-11.0,FEELPP_CXX=icpc
#
####################################################################

# process the arguments

set(ARGLIST ${CTEST_SCRIPT_ARG})
while(${ARGLIST} MATCHES  ".+.*")

  # pick first
  string(REGEX MATCH "([^,]*)(,.*)?" DUMMY ${ARGLIST})
  SET(TOP ${CMAKE_MATCH_1})

  # remove first
  string(REGEX MATCHALL "[^,]*,(.*)" DUMMY ${ARGLIST})
  SET(ARGLIST ${CMAKE_MATCH_1})

  # decompose as a pair key=value
  string(REGEX MATCH "([^=]*)(=.*)?" DUMMY ${TOP})
  SET(KEY ${CMAKE_MATCH_1})

  string(REGEX MATCH "[^=]*=(.*)" DUMMY ${TOP})
  SET(VALUE ${CMAKE_MATCH_1})

  # set the variable to the specified value
  if(VALUE)
    SET(${KEY} ${VALUE})
  else(VALUE)
    SET(${KEY} ON)
  endif(VALUE)

endwhile(${ARGLIST} MATCHES ".+.*")

####################################################################
# Automatically set some user variables if they have not been defined manually
####################################################################
cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD}  ARGS -s OUTPUT_VARIABLE SITE_HOSTNAME)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

#site_name calls hostname and the the result in FEELPP_SITE
if(NOT FEELPP_SITE)
  site_name(FEELPP_SITE)
endif(NOT FEELPP_SITE)

#Read the configuration for cdash, update process (svn/git) and subprojects to build.
if ( EXISTS ${FEELPP_CTEST_CONFIG} )
  include(${FEELPP_CTEST_CONFIG})
  message("FEELPP_CTEST_CONFIG: " ${FEELPP_CTEST_CONFIG})
  set(FEELPP_BUILD_STRING "${OS_VERSION}-${ARCH}")
endif()

#Check the compiler
if (${FEELPP_CXXNAME} MATCHES "gcc*")
  message("GCC")
  if (DEFINED GCC_MAKE_ARGS)
    set(MAKE_ARGS "${GCC_MAKE_ARGS}")
  endif()
  if (DEFINED GCC_PARALLEL)
    set(PARALLEL "${GCC_PARALLEL}")
  endif()
elseif (${FEELPP_CXXNAME} MATCHES "clang")
  message("clang")
  if (DEFINED CLANG_MAKE_ARGS)
    set(MAKE_ARGS "${CLANG_MAKE_ARGS}")
  endif()
  if (DEFINED CLANG_PARALLEL)
    set(PARALLEL "${CLANG_PARALLEL}")
  endif()
endif()
set(CTEST_BUILD_FLAGS -j${PARALLEL})
set(CTEST_PARALLEL_LEVEL ${PARALLEL})
message("MAKE_ARGS -- ${MAKE_ARGS}")
message("PARALLEL -- ${PARALLEL}")


if(NOT FEELPP_CMAKE_DIR)
  SET(FEELPP_CMAKE_DIR "")
endif(NOT FEELPP_CMAKE_DIR)

if (NOT FEELPP_CXX)
  MESSAGE(WARNING "NO CXX COMPILER PROVIDED - USING DEFAULT ONE")
  set(FEELPP_CXX "g++")
endif(NOT FEELPP_CXX)

#Generate default build string OS-VERSION-ARCH-COMPILER
if(NOT FEELPP_BUILD_STRING)
  # OS
  build_name(FEELPP_OS_VERSION)

  # arch
  set(FEELPP_ARCH ${CMAKE_SYSTEM_PROCESSOR})
  if(WIN32)
    set(FEELPP_ARCH $ENV{PROCESSOR_ARCHITECTURE})
  else(WIN32)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE FEELPP_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif(WIN32)

  set(FEELPP_BUILD_STRING ${FEELPP_OS_VERSION}${FEELPP_ARCH}-${FEELPP_CXXNAME})

endif(NOT FEELPP_BUILD_STRING)

if(DEFINED FEELPP_EXPLICIT_VECTORIZATION)
  set(FEELPP_BUILD_STRING ${FEELPP_BUILD_STRING}-${FEELPP_EXPLICIT_VECTORIZATION})
endif(DEFINED FEELPP_EXPLICIT_VECTORIZATION)

if(NOT FEELPP_WORK_DIR)
  set(FEELPP_WORK_DIR ${CTEST_SCRIPT_DIRECTORY})
endif(NOT FEELPP_WORK_DIR)

if(NOT CTEST_SOURCE_DIRECTORY)
  SET (CTEST_SOURCE_DIRECTORY "${FEELPP_WORK_DIR}/feelpp")
endif(NOT CTEST_SOURCE_DIRECTORY)

if(NOT FEELPP_MODE)
  set(FEELPP_MODE Nightly)
endif(NOT FEELPP_MODE)

if(NOT CTEST_BINARY_DIRECTORY)
  SET (CTEST_BINARY_DIRECTORY "${FEELPP_WORK_DIR}/${FEELPP_MODE}_${FEELPP_CXXNAME}")
endif(NOT CTEST_BINARY_DIRECTORY)
## Remove ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt - usefull if config has changed (petsc, boost...)
execute_process(COMMAND "rm ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
message(STATUS "CMakeCache.txt removed")

## mandatory variables (the default should be ok in most cases):

#if(NOT FEELPP_NO_UPDATE)
find_program(CTEST_GIT_COMMAND NAMES git)
#find_program(CTEST_SVN_COMMAND NAMES svn)
#SET (CTEST_SVN_CHECKOUT   "${CTEST_GIT_COMMAND} co svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk ${CTEST_SOURCE_DIRECTORY}")
#SET (CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://code.google.com/p/feelpp/")
#SET (CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone https://github.com/feelpp/feelpp.git")
set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

#########################################
# UPDATE source dir & modules if defined
#########################################
#execute_process(COMMAND "cd ${CTEST_SOURCE_DIRECTORY} && ${CTEST_GIT_COMMAND} pull" OUTPUT_VARIABLE MODULE_OUTPUT)
#message(STATUS "updated ${CTEST_SOURCE_DIRECTORY} : ${MODULE_OUTPUT}")
foreach(module ${FEELPP_MODULES})
  # update the modules using svn update
  execute_process(
    #COMMAND "cd ${CTEST_SOURCE_DIRECTORY}/${module} && ${CTEST_SVN_COMMAND} update"
    COMMAND "cd ${CTEST_SOURCE_DIRECTORY}/${module} && ${CTEST_GIT_COMMAND} pull"
    OUTPUT_VARIABLE MODULE_OUTPUT)
  message(STATUS "updated ${module} : ${MODULE_OUTPUT}")
endforeach()
####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this make sure we get consistent outputs
SET($ENV{LC_MESSAGES} "en_EN")

if (UNIX)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif(UNIX)

#The idea behind ctest launchers is that they wrap each compile or link step so
#the output can be saved and sent to CDash in the event of a warning or
#error. Rather than trying to grep through and analyze the full build output
#after thousands of compile and link calls, with this technique, ctest may
#simply capture the error output directly and pass it in its entirety to the
#dashboard. This helps immensely in figuring out some why some errors occur,
#without necessarily even having access to the client machine.
set(CTEST_USE_LAUNCHERS 1)
find_program(CTEST_CMAKE_COMMAND NAMES "cmake")
message("CMAKE FOUND -- ${CTEST_CMAKE_COMMAND}")
# Generating the CTEST_CONFIGURE_COMMAND variable
set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND}")
set(CMAKE_OPTIONS "-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS}")
set(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DCMAKE_CXX_COMPILER:STRING=${FEELPP_CXX} -DCMAKE_C_COMPILER:STRING=/usr/bin/gcc")
set(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DVTK_LEGACY_REMOVE:BOOL=ON")

##
#please see: http://public.kitware.com/pipermail/cdash/2010-June/000820.html
#note: CTEST_CMAKE_COMMAND is intended to be a read-only variable
##
if(DEFINED FEELPP_CMAKE_ARGS)
  set(CMAKE_OPTIONS "${CMAKE_OPTIONS} ${FEELPP_CMAKE_ARGS}")
endif(DEFINED FEELPP_CMAKE_ARGS)
#Define the cmake options
if(DEFINED BUILD_TYPE )
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DCMAKE_BUILD_TYPE=${BUILD_TYPE}")
endif()
if(DEFINED ENABLE_ALTIVEC)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_ALTIVEC=${ENABLE_ALTIVEC}")
endif()
if(DEFINED ENABLE_BUILD_STATIC)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_BUILD_STATIC=${ENABLE_BUILD_STATIC}")
endif()
if(DEFINED ENABLE_DOXYGEN)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_DOXYGEN=${ENABLE_DOXYGEN}")
endif()
if(DEFINED ENABLE_NEON)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_NEON=${ENABLE_NEON}")
endif()
if(DEFINED ENABLE_OPENTURNS)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_OPENTURNS=${ENABLE_OPENTURNS}")
endif()
if(DEFINED ENABLE_PCH_FOR_APPLICATIONS)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_PCH_FOR_APPLICATIONS=${ENABLE_PCH_FOR_APPLICATIONS}")
endif()
if(DEFINED ENABLE_VERBOSE_CMAKE)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_VERBOSE_CMAKE=${ENABLE_VERBOSE_CMAKE}")
endif()
if(DEFINED BENCHMARK_FLAG)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_BENCHMARK_FLAG=${BENCHMARK_FLAG}")
endif()
if(DEFINED ENABLE_TESTS)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_TESTS=${ENABLE_TESTS}")
endif()
if(DEFINED ENABLE_DOCUMENTATION)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_DOCUMENTATION=${ENABLE_DOCUMENTATION}")
endif()
if(DEFINED ENABLE_BENCHMARKS)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_BENCHMARKS=${ENABLE_BENCHMARKS}")
endif()
if(DEFINED ENABLE_RESEARCH)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_RESEARCH=${ENABLE_RESEARCH}")
endif()
if(DEFINED ENABLE_APPLICATIONS)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_APPLICATIONS=${ENABLE_APPLICATIONS}")
endif()
if(DEFINED ENABLE_CRB_ALL)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_CRB_ALL=${ENABLE_CRB_ALL}")
endif()
if(DEFINED ENABLE_APPLICATIONS_CRB)
  SET(CMAKE_OPTIONS "${CMAKE_OPTIONS} -DFEELPP_ENABLE_APPLICATIONS_CRB=${ENABLE_APPLICATIONS_CRB}")
endif()

# set test timeout to 300s
set(CTEST_TIMEOUT "300")

# raise the warning/error limit
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "33331")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "33331")

# to get CTEST_PROJECT_SUBPROJECTS definition:
include("${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake")

# clear the binary directory and create an initial cache
#CTEST_EMPTY_BINARY_DIRECTORY (${CTEST_BINARY_DIRECTORY})
set(CTEST_INITIAL_CACHE "CMAKE_CXX_COMPILER:STRING=${FEELPP_CXX}")
# site
set(CTEST_SITE "${FEELPP_SITE}")

set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/cmake/modules)
# build name
set(_BLOCK_ true)
FIND_PACKAGE(Boost)
FIND_PACKAGE(PETSc)
set(CTEST_BUILD_NAME "${FEELPP_BUILD_STRING}-${FEELPP_CXXNAME}-petsc_${PETSC_VERSION}-boost_${Boost_VERSION}")
# string(REPLACE "+" "%2B" CTEST_BUILD_NAME ${CTEST_BUILD_NAME})
# set(CTEST_BUILD_NAME "${FEELPP_BUILD_STRING}-${FEELPP_CXXNAME}")
# should ctest wipe the binary tree before running
#SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

if(FEELPP_CXX AND NOT WIN32)
  set(CTEST_ENVIRONMENT "CXX=${FEELPP_CXX}")
endif(FEELPP_CXX AND NOT WIN32)

ctest_start(${FEELPP_MODE})
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}") #git pull

## CMAKE has to be called once !
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CMAKE_OPTIONS} ${CTEST_SOURCE_DIRECTORY}")
ctest_configure() # Execute CTEST_CONFIGURE_COMMAND
# Submit results to a dashboard server.
ctest_submit(PARTS Configure)
ctest_submit(PARTS Update Notes)

foreach(subproject ${CTEST_PROJECT_SUBPROJECTS})
  set_property(GLOBAL PROPERTY SubProject ${subproject})
  set_property(GLOBAL PROPERTY Label ${subproject})

  message(WARNING "build target " ${subproject})
  #set(CTEST_BUILD_TARGET “${subproject}”)
  set(CTEST_BUILD_TARGET "${subproject}")
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
  ctest_submit(PARTS Build)

  message(WARNING "BUILD "${CTEST_BINARY_DIRECTORY})
  # runs only tests that have a LABELS property matching "${subproject}"
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" INCLUDE_LABEL "${subproject}"  )
  # Submit results to a dashboard server.
  ctest_submit(PARTS Test)
endforeach()
