cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

project(CLN)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/cln/version.h" _cl_vinfo REGEX "^#define[\t ]+CL_VERSION_.*")
string(REGEX REPLACE "^.*CL_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" CL_VERSION_MAJOR "${_cl_vinfo}")
string(REGEX REPLACE "^.*CL_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" CL_VERSION_MINOR "${_cl_vinfo}")
string(REGEX REPLACE "^.*CL_VERSION_PATCHLEVEL[ \t]+([0-9]+).*" "\\1" CL_VERSION_PATCHLEVEL "${_cl_vinfo}")
set(CL_VERSION "${CL_VERSION_MAJOR}.${CL_VERSION_MINOR}.${CL_VERSION_PATCHLEVEL}")

# Library
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/include/cln/version.h _cl_vinfo REGEX "^#define[\t ]+CL_LT_.*")
string(REGEX REPLACE "^.*CL_LT_CURRENT[ \t]+([0-9]+).*" "\\1" cl_lt_current "${_cl_vinfo}")
string(REGEX REPLACE "^.*CL_LT_AGE[ \t]+([0-9]+).*" "\\1" cl_lt_age "${_cl_vinfo}")
string(REGEX REPLACE "^.*CL_LT_REVISION[ \t]+([0-9]+).*" "\\1" cl_lt_revision "${_cl_vinfo}")

# XXX: CMake has no portable library versioning?
math(EXPR libcln_soversion "${cl_lt_current} - ${cl_lt_age}")
set(libcln_version ${libcln_soversion}.${cl_lt_age}.${cl_lt_revision})
option(CLN_USE_GMP "Use GMP low level routines" ON)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config.status
		OR EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/libtool)
	message(FATAL_ERROR
		"Source directory has been already configured by configure script. "
		"Please run `make distclean` in ${CMAKE_CURRENT_SOURCE_DIR} and try again."
	)
endif()
include(GNUInstallDirs)

set(cl_config_files
	include/cln/config.h
	include/cln/host_cpu.h
	include/cln/intparam.h
	autoconf/cl_config.h
	src/base/cl_base_config.h
	src/base/cl_gmpconfig.h
	src/timing/cl_t_config.h
)
set(cln_generated_headers)

include(CheckIncludeFile)
include(CheckLibraryExists)
include(CheckTypeSize)
include(TestBigEndian)
include(CheckFunctionExists)
include(CheckCXXSourceCompiles)
include(CheckTypeAlign)

enable_testing()
if (NOT TARGET check)
	add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()
if (NOT TARGET test_suite)
	add_custom_target(test_suite)
endif()
if (NOT TARGET info)
	add_custom_target(info ALL)
endif()
if (NOT TARGET html)
	add_custom_target(html)
endif()
if (NOT TARGET pdf)
	add_custom_target(pdf)
endif()

if (NOT DEFINED BUILD_SHARED_LIBS)
	if (NOT MSVC)
		set(BUILD_SHARED_LIBS true)
	else()
		set(BUILD_SHARED_LIBS false)
	endif()
endif()
if (WIN32)
	if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
		set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
	endif()
endif()

if (CLN_USE_GMP)
find_package(GMP REQUIRED)
endif()

if (GMP_FOUND)
	set(CL_USE_GMP TRUE)
endif()

if (GMP_FOUND)
	get_filename_component(GMP_LIBDIR ${GMP_LIBRARIES} PATH)
	include(GMPLimbSize)
endif(GMP_FOUND)


check_type_size("short" SIZEOF_SHORT)
check_type_size("int" SIZEOF_INT)
check_type_size("long" SIZEOF_LONG)
check_type_size("long long" SIZEOF_LONGLONG)
check_type_size("void *" SIZEOF_VOIDP)
check_type_align("void *" ALIGNOF_VOIDP)

# not guaranteed by the C standard, but holds in practice
set(cl_char_bitsize 8)
math(EXPR x "${SIZEOF_SHORT} * ${cl_char_bitsize}")
set(cl_short_bitsize ${x} CACHE INTERNAL "short_bitsize")
math(EXPR cl_int_bitsize "${SIZEOF_INT} * ${cl_char_bitsize}")
math(EXPR cl_long_bitsize "${SIZEOF_LONG} * ${cl_char_bitsize}")
math(EXPR cl_long_long_bitsize "${SIZEOF_LONGLONG} * ${cl_char_bitsize}")
math(EXPR cl_pointer_bitsize "${SIZEOF_VOIDP} * ${cl_char_bitsize}")


check_include_file("unistd.h" HAVE_UNISTD_H)
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)

set(_save_required_includes ${CMAKE_REQUIRED_INCLUDES})
set(_save_required_libraries ${CMAKE_REQUIRED_LIBRARIES})
set(_save_required_flags ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-Werror")
set(CMAKE_REQUIRED_INCLUDES "")
set(CMAKE_REQUIRED_LIBRARIES "")
CHECK_CXX_SOURCE_COMPILES("
	void f() __attribute__((flatten));
	int main() { return 0; }
	"
	CL_HAVE_ATTRIBUTE_FLATTEN
)
set(CMAKE_REQUIRED_FLAGS ${_save_required_flags})
set(CMAKE_REQUIRED_INCLUDES ${_save_required_includes})
set(CMAKE_REQUIRED_LIBRARIES ${_save_required_libraries})


test_big_endian(cl_cv_bigendian_p)
if(cl_cv_bigendian_p)
	set(short_big_endian 1 CACHE INTERNAL "CPU is big endian")
	set(int_big_endian 1 CACHE INTERNAL "CPU is big endian")
	set(long_big_endian 1 CACHE INTERNAL "CPU is big endian")
	set(long_long_big_endian 1 CACHE INTERNAL "CPU is big endian")
else()
	set(short_little_endian 1 CACHE INTERNAL "CPU is little endian")
	set(int_little_endian 1 CACHE INTERNAL "CPU is little endian")
	set(long_little_endian 1 CACHE INTERNAL "CPU is little endian")
	set(long_long_little_endian 1 CACHE INTERNAL "CPU is little endian")
endif()

include(AsmUnderscore)
include(AsmNoexecstack)

set(GMP_LIBDIR_PC)

set(_cln_rpath_reloc "$ORIGIN/../lib")
set(_wl_rpath "${CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG}")
set(_wl_rpath_link "${CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG}")

# rpath for the pkg-config meta-data
set(_libcln_rpath "${_cln_rpath_reloc}")
set(_libcln_pc_rpath "${_wl_rpath}${_cln_rpath_reloc}")

list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
	set(_libcln_pc_rpath "${_libcln_pc_rpath} ${_wl_rpath}\${libdir}")
	list(APPEND _libcln_rpath ${CMAKE_INSTALL_FULL_LIBDIR})
endif()

list(FIND CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${GMP_LIBDIR}" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
	list(APPEND _libcln_rpath ${GMP_LIBDIR})
	set(_libcln_pc_rpath "${_libcln_pc_rpath} ${_wl_rpath_link}${GMP_LIBDIR}")
	set(GMP_LIBDIR_PC "-L${GMP_LIBDIR}")
endif()
string(REPLACE ";" ":" libcln_rpath "${_libcln_rpath}")

if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_RPATH)
	set(CMAKE_INSTALL_RPATH_USE_LINK_RPATH TRUE)
endif()
if (NOT DEFINED CMAKE_INSTALL_RPATH)
	set(CMAKE_INSTALL_RPATH ${libcln_rpath})
	if (_wl_rpath)
		set(CLN_PC_RPATH ${_libcln_pc_rpath})
	endif()
endif()
if (APPLE AND NOT DEFINED CMAKE_INSTALL_NAME_DIR)
	set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR})
endif()

macro(cl_config_file relname)
	set(_dst ${CMAKE_CURRENT_BINARY_DIR}/${relname})
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${relname}.cmake ${_dst})
	if ("${relname}" MATCHES "^include")
		list(APPEND cln_generated_headers "${_dst}")
	endif()
endmacro()

foreach(_cf ${cl_config_files})
	cl_config_file("${_cf}")
endforeach()
install(FILES ${cln_generated_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cln")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cln.pc.cmake ${CMAKE_CURRENT_BINARY_DIR}/cln.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cln.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/cln-config-version.cmake
	VERSION ${CL_VERSION}
	COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/cln-config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/cln-config.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cln
)
install(FILES
	${CMAKE_CURRENT_LIST_DIR}/cmake/modules/FindGMP.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cln-config.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cln-config-version.cmake
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cln
)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(include)
add_subdirectory(examples)
add_subdirectory(doc)

