cmake_minimum_required(VERSION 3.10)
project(CLN_examples)

set(cln_examples
	contfrac
	e
	fibonacci
	legendre
	lucaslehmer
	nextprime
	perfnum
	pi
)

if (NOT TARGET cln)
	find_package(cln REQUIRED CONFIG)
endif()
macro(cl_add_example srcname)
	add_executable(${srcname} ${srcname}.cc)
	target_link_libraries(${srcname} cln::cln)
endmacro()

foreach(_ex ${cln_examples})
	cl_add_example(${_ex})
endforeach()

install(TARGETS pi RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

