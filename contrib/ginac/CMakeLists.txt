cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

project(GiNaC)
file(STRINGS ${CMAKE_SOURCE_DIR}/contrib/ginac/ginac/version.h _ginac_vinfo REGEX "^#define[\t ]+GINACLIB_.*_VERSION.*")
string(REGEX REPLACE "^.*GINACLIB_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" GINAC_MAJOR_VERSION "${_ginac_vinfo}")
string(REGEX REPLACE "^.*GINACLIB_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" GINAC_MINOR_VERSION "${_ginac_vinfo}")
string(REGEX REPLACE "^.*GINACLIB_MICRO_VERSION[ \t]+([0-9]+).*" "\\1" GINAC_MICRO_VERSION "${_ginac_vinfo}")
set(GINAC_VERSION "${GINAC_MAJOR_VERSION}.${GINAC_MINOR_VERSION}.${GINAC_MICRO_VERSION}")

# Library versioning info
file(STRINGS ${CMAKE_SOURCE_DIR}/contrib/ginac/ginac/version.h _ginac_vinfo REGEX "^#define[\t ]+GINAC_LT_.*")
string(REGEX REPLACE "^.*GINAC_LT_CURRENT[ \t]+([0-9]+).*" "\\1" ginac_lt_current "${_ginac_vinfo}")
string(REGEX REPLACE "^.*GINAC_LT_AGE[ \t]+([0-9]+).*" "\\1" ginac_lt_age "${_ginac_vinfo}")
string(REGEX REPLACE "^.*GINAC_LT_REVISION[ \t]+([0-9]+).*" "\\1" ginac_lt_revision "${_ginac_vinfo}")
# XXX: it looks like we need to set this for every platform?
math(EXPR ginaclib_soversion "${ginac_lt_current} - ${ginac_lt_age}")
set(ginaclib_version ${ginaclib_soversion}.${ginac_lt_age}.${ginac_lt_revision})

# # for supermuc this option should be OFF
# option(USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER "Enable use of standard include headers in ginac-excompiler" ON)
IF(USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER)
  #add -DUSE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER to CFLAGS
  message( STATUS "[contrib ginac] Enable use of standard headers")
  set(_USE_STANDARD_C_HEADERS_ "-DUSE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER")
ENDIF()

# make check
enable_testing()
add_custom_target(check-ginac COMMAND ${CMAKE_CTEST_COMMAND})
# make info
#add_custom_target(info ALL)
#add_custom_target(html)
#add_custom_target(pdf)

find_package(CLN 1.3.2)
#if ( CLN_FOUND )
include_directories(${CLN_INCLUDE_DIR})
#endif()

include(CheckIncludeFile)
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("unistd.h" HAVE_UNISTD_H)

# Check Python version
# If python is version 3, we need to patch the scripts/yaptu.py file
# (exec is a statement in Python 2 and a function in Python 3)
FIND_PACKAGE(PythonInterp)

if(PYTHONINTERP_FOUND)
    if(PYTHON_VERSION_MAJOR MATCHES "3")
        message(STATUS "Python 3 detected: patching files")
        get_filename_component(PYTHONINTERP_PATH ${PYTHON_EXECUTABLE} PATH)
        execute_process(COMMAND ${PYTHONINTERP_PATH}/2to3 -w --no-diffs ${CMAKE_SOURCE_DIR}/contrib/ginac/)
        message(WARNING "Ginac files have been patched for Python 3 compatibility")
    elseif(PYTHON_VERSION_MAJOR MATCHES "2")
        message(STATUS "Python 2 detected: No need to patch")
    else()
        message(FATAL_ERROR "Unknown Python version: ${PYTHON_VERSION_STRING}")
    endif()
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/ginac)

# This macro implements some very special logic how to deal with the cache.
# By default the various install locations inherit their value from their
#"parent" variable, so if you set CMAKE_INSTALL_PREFIX, then
# EXEC_INSTALL_PREFIX, BIN_INSTALL_DIR, LIB_INSTALL_DIR, etc will calculate
# their value by appending subdirs to CMAKE_INSTALL_PREFIX.
# This would work just fine without using the cache.
# But if somebody wants e.g. a different EXEC_INSTALL_PREFIX this value
# has to go into the cache, otherwise it will be forgotten on the next cmake
# run. Once a variable is in the cache, it doesn't depend on its "parent"
# variables anymore and you can only change it by editing it directly.
# This macro helps in this regard, because as long as you don't set one
# of the variables explicitly to some location, the value will be computed
# from parents of the variable in question. So modifying CMAKE_INSTALL_PREFIX
# later on will have the desired effect.
# But once you decide to set e.g. EXEC_INSTALL_PREFIX to some special
# location this will go into the cache and it will no longer depend on
# CMAKE_INSTALL_PREFIX.

macro(_set_fancy _var _value _comment)
	set(predefinedvalue "${_value}")
	if ("${CMAKE_INSTALL_PREFIX}" STREQUAL "${GINAC_INSTALL_DIR}" AND DEFINED GINAC_${_var})
		set(predefinedvalue "${GINAC_${_var}}")
	endif()
	if (NOT DEFINED ${_var})
		set(${_var} ${predefinedvalue})
	else()
		set(${_var} "${${_var}}" CACHE PATH "${_comment}")
	endif()
endmacro(_set_fancy)

_set_fancy(EXEC_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}"
	   "Base directory for libraries and executables")
_set_fancy(LIB_INSTALL_DIR "${EXEC_INSTALL_PREFIX}/lib"
	   "Libraries installation directory")
_set_fancy(BIN_INSTALL_DIR "${EXEC_INSTALL_PREFIX}/bin"
	   "Binaries installation directory")
_set_fancy(SHARE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/share"
	   "Base directory for architecture independent files")
_set_fancy(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include"
	   "Headers installation directory")

if (NOT DEFINED BUILD_SHARED_LIBS)
	if (NOT MSVC)
		set(BUILD_SHARED_LIBS true)
	else()
		set(BUILD_SHARED_LIBS false)
	endif()
endif()

# Set proper rpath so tools will Just Work(TM) after make install.
# Also take care to add -Wl,-rpath, stanza into the *.pc file so that
#
# g++ `pkg-config --cflags --libs ginac`
#
# will Just Work (TM), too.
# Distro packagers should use -DCMAKE_INSTALL_RPATH="" to avoid
# setting rpath on installed binaries.

# rpath for making binaries/libraries relocatable
set(_ginac_rpath_reloc "$ORIGIN/../lib")
set(_wl_rpath "${CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG}")

# rpath for the pkg-config meta-data.
set(_ginaclib_rpath "${_wl_rpath}${_ginac_rpath_reloc}")
list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${LIB_INSTALL_DIR}" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
	list(APPEND _ginaclib_rpath "${_wl_rpath}\${libdir}")
endif()
set(GINACLIB_RPATH)
if (NOT CMAKE_SKIP_RPATH)
	if (_wl_rpath)
		set(GINACLIB_RPATH "${_ginaclib_rpath}")
	endif()
endif()

configure_file(${CMAKE_SOURCE_DIR}/contrib/ginac/ginac.pc.cmake ${CMAKE_BINARY_DIR}/ginac.pc @ONLY)
#install(FILES ${CMAKE_BINARY_DIR}/ginac.pc DESTINATION "${LIB_INSTALL_DIR}/pkgconfig")

# rpath for libginac.so itself, ginsh, and friends
set(_ginac_rpath ${_ginac_rpath_reloc})
foreach(_d ${CLN_LIBRARY_DIRS} ${LIB_INSTALL_DIR})
	list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${_d}" isSystemDir)
	if ("${isSystemDir}" STREQUAL "-1")
		list(APPEND _ginac_rpath "${_d}")
	endif()
endforeach()
list(REMOVE_DUPLICATES _ginac_rpath)
string(REPLACE ";" ":" ginac_rpath "${_ginac_rpath}")

if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_RPATH)
	set(CMAKE_INSTALL_RPATH_USE_LINK_RPATH TRUE)
endif()
if (NOT DEFINED CMAKE_INSTALL_RPATH)
	set(CMAKE_INSTALL_RPATH ${ginac_rpath})
endif()
if (APPLE AND NOT DEFINED CMAKE_INSTALL_NAME_DIR)
	set(CMAKE_INSTALL_NAME_DIR ${LIB_INSTALL_DIR})
endif()

list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${LIB_INSTALL_DIR}" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
	string(REPLACE ":" ";" _install_rpath "${CMAKE_INSTALL_RPATH}")
	list(FIND _install_rpath "${LIB_INSTALL_DIR}" _is_rpath_consistent)
	if ("${_is_rpath_consistent}" STREQUAL "-1")
		message(WARNING "the libginac.so library will be installed into "
				"a non-standard directory (${LIB_INSTALL_DIR}), "
				"however, the rpath (${_install_rpath}) "
				"does not contain that directory. Most likely "
				"things won't work without extra configuration "
				"(tweaking LD_LIBRARY_PATH, /etc/ld.so.conf, etc).")
	endif()
endif()


include(FindFLEX)
include(FindBISON)
find_package(BISON)
find_package(FLEX)
find_package(Readline)
if (READLINE_FOUND)
	set(HAVE_LIBREADLINE 1)
	set(HAVE_READLINE_READLINE_H 1)
	set(HAVE_READLINE_HISTORY_H 1)
endif()

find_program(MAKEINFO makeinfo)
find_program(FIG2DEV fig2dev)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(ginac)
add_subdirectory(tools)
# feel++ edit
#add_subdirectory(check)
if (BISON_FOUND AND FLEX_FOUND)
	add_subdirectory(ginsh)
endif()
# feel++ edit
#if (MAKEINFO)
#	add_subdirectory(doc)
#endif()
