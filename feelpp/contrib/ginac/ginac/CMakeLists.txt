cmake_minimum_required(VERSION 3.13)

configure_file( excompiler.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/excompiler.cpp)

set(ginaclib_sources
    add.cpp
    archive.cpp
    basic.cpp
    clifford.cpp
    color.cpp
    constant.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/excompiler.cpp
    ex.cpp
    expair.cpp
    expairseq.cpp
    exprseq.cpp
    factor.cpp
    fail.cpp
    fderivative.cpp
    function.cpp
    idx.cpp
    indexed.cpp
    inifcns.cpp
    inifcns_gamma.cpp
    inifcns_nstdsums.cpp
    inifcns_trans.cpp
    integral.cpp
    lst.cpp
    matrix.cpp
    mul.cpp
    ncmul.cpp
    normal.cpp
    numeric.cpp
    operators.cpp
    parser/default_reader.cpp
    parser/lexer.cpp
    parser/parse_binop_rhs.cpp
    parser/parse_context.cpp
    parser/parser_compat.cpp
    parser/parser.cpp
    polynomial/chinrem_gcd.cpp
    polynomial/collect_vargs.cpp
    polynomial/cra_garner.cpp
    polynomial/divide_in_z_p.cpp
    polynomial/gcd_uvar.cpp
    polynomial/mgcd.cpp
    polynomial/mod_gcd.cpp
    polynomial/normalize.cpp
    polynomial/optimal_vars_finder.cpp
    polynomial/pgcd.cpp
    polynomial/primpart_content.cpp
    polynomial/remainder.cpp
    polynomial/upoly_io.cpp
    power.cpp
    print.cpp
    pseries.cpp
    registrar.cpp
    relational.cpp
    remember.cpp
    symbol.cpp
    symmetry.cpp
    tensor.cpp
    utils.cpp
    wildcard.cpp
)

set(ginaclib_public_headers
    ginac.h
    add.h
    archive.h
    assertion.h
    basic.h
    class_info.h
    clifford.h
    color.h
    compiler.h
    constant.h
    container.h
    ex.h
    excompiler.h
    expair.h
    expairseq.h
    exprseq.h
    fail.h
    factor.h
    fderivative.h
    flags.h
    ${CMAKE_CURRENT_BINARY_DIR}/function.h
    hash_map.h
    idx.h
    indexed.h
    inifcns.h
    integral.h
    lst.h
    matrix.h
    mul.h
    ncmul.h
    normal.h
    numeric.h
    operators.h
    power.h
    print.h
    pseries.h
    ptr.h
    registrar.h
    relational.h
    structure.h
    symbol.h
    symmetry.h
    tensor.h
    version.h
    wildcard.h
    parser/parser.h
    parser/parse_context.h
)

set(ginaclib_private_headers
    remember.h
    utils.h
    crc32.h
    hash_seed.h
    parser/lexer.h
    parser/debug.h
    polynomial/gcd_euclid.h
    polynomial/remainder.h
    polynomial/normalize.h
    polynomial/upoly.h
    polynomial/ring_traits.h
    polynomial/mod_gcd.h
    polynomial/cra_garner.h
    polynomial/upoly_io.h
    polynomial/prem_uvar.h
    polynomial/eval_uvar.h
    polynomial/interpolate_padic_uvar.h
    polynomial/sr_gcd_uvar.h
    polynomial/heur_gcd_uvar.h
    polynomial/chinrem_gcd.h
    polynomial/collect_vargs.h
    polynomial/divide_in_z_p.h
    polynomial/euclid_gcd_wrap.h
    polynomial/eval_point_finder.h
    polynomial/newton_interpolate.h
    polynomial/optimal_vars_finder.h
    polynomial/pgcd.h
    polynomial/poly_cra.h
    polynomial/primes_factory.h
    polynomial/smod_helpers.h
    polynomial/debug.h
)

add_library(feelpp_ginac ${ginaclib_sources})
set_property(TARGET feelpp_ginac PROPERTY CXX_STANDARD ${FEELPP_STD_CPP})
add_definitions(-DLIBEXECDIR="${LIBEXECDIR}/")
set_target_properties(feelpp_ginac PROPERTIES
	SOVERSION ${ginaclib_soversion}
	VERSION ${ginaclib_version})
target_link_libraries(feelpp_ginac cln::cln ${CMAKE_DL_LIBS} ${Boost_MPI_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${MPI_LIBRARIES})
target_include_directories(feelpp_ginac
  PUBLIC
  # ${CLN_INCLUDE_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/feelpp/>
  PRIVATE
  ${MPI_CXX_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS})

if (NOT BUILD_SHARED_LIBS)
    set_target_properties(feelpp_ginac PROPERTIES OUTPUT_NAME "feelpp_ginac")
	# Avoid the static library vs import library conflict (necessary for
	# m$ toolchain).
	set_target_properties(feelpp_ginac PROPERTIES PREFIX "lib")
  else()
    set_property(TARGET feelpp_ginac PROPERTY MACOSX_RPATH ON)
endif()

install(TARGETS feelpp_ginac EXPORT feelpp-contrib-export-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Libs
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Libs
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Libs)
install(FILES ${ginaclib_public_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ginac" COMPONENT Devel)

## Feel++ modification
## Those headers will be installed twice, Needed for custom installations to properly work and not breaking codes for now
set(ginaclib_public_parser_headers
      parser/parser.h
      parser/parse_context.h)
install(FILES ${ginaclib_public_parser_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ginac/parser" COMPONENT Devel)

# # End Feel++ modification
if( Python3_EXECUTABLE )
  set(GINAC_PYTHON_EXECUTABLE ${Python3_EXECUTABLE} )
else()
  find_program(GINAC_PYTHON_EXECUTABLE NAMES python)
  if ( NOT GINAC_PYTHON_EXECUTABLE )
    message(FATAL_ERROR "[Ginac] not found python executable)")
  endif()
endif()

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/function.h
	COMMAND ${GINAC_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/function.py -o ${CMAKE_CURRENT_BINARY_DIR}/function.h ${CMAKE_CURRENT_SOURCE_DIR}/function.hppy
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/function.hppy ${CMAKE_CURRENT_SOURCE_DIR}/function.py)

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/function.cpp
	COMMAND ${GINAC_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/function.py -o ${CMAKE_CURRENT_BINARY_DIR}/function.cpp ${CMAKE_CURRENT_SOURCE_DIR}/function.cppy
	DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/function.h ${CMAKE_CURRENT_SOURCE_DIR}/function.cppy ${CMAKE_CURRENT_SOURCE_DIR}/function.py)
