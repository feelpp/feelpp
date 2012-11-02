include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ginac)
add_definitions(-DIN_GINAC)

set(ginac_tests
	check_numeric
	check_inifcns
	check_matrices
	check_lsolve
	heur_gcd_bug
	exam_paranoia
	exam_heur_gcd
	match_bug
	parser_bugs
	exam_numeric_archive
	exam_numeric
	exam_powerlaws
	exam_inifcns
	exam_inifcns_nstdsums
	exam_differentiation
	exam_polygcd
	exam_normalization
	exam_factor
	exam_pseries
	exam_matrices
	exam_lsolve
	exam_indexed
	exam_color
	exam_clifford
	exam_archive
	exam_structure
	exam_hashmap
	exam_misc
	exam_mod_gcd
	exam_cra
	bugme_chinrem_gcd
	pgcd_relatively_prime_bug
	pgcd_infinite_loop)

set(ginac_timings
	time_dennyfliegner
	time_gammaseries
	time_vandermonde
	time_toeplitz
	time_hashmap
	time_lw_A
	time_lw_B
	time_lw_C
	time_lw_D
	time_lw_E
	time_lw_F
	time_lw_G
	time_lw_H
	time_lw_IJKL
	time_lw_M1
	time_lw_M2
	time_lw_N
	time_lw_O
	time_lw_P
	time_lw_Pprime
	time_lw_Q
	time_lw_Qprime
	time_antipode
	time_fateman_expand
	time_uvar_gcd
	time_parser)

macro(add_ginac_test thename)
	if ("${${thename}_sources}" STREQUAL "")
		set(${thename}_sources ${thename}.cpp ${${thename}_extra_src})
	endif()
	add_executable(${thename} EXCLUDE_FROM_ALL ${${thename}_sources})
	target_link_libraries(${thename} ginac)
	add_dependencies(check-ginac ${thename})
	add_test(NAME ${thename} COMMAND ${thename}${CMAKE_EXECUTABLE_SUFFIX})
endmacro()

macro(add_ginac_timing thename)
	set(${thename}_extra_src timer.cpp randomize_serials.cpp)
	add_ginac_test(${thename})
endmacro()

set(check_matrices_extra_src genex.cpp)
set(check_lsolve_extra_src genex.cpp)
set(exam_heur_gcd_sources heur_gcd_bug.cpp)
set(exam_numeric_archive_sources numeric_archive.cpp)

foreach(tst ${ginac_tests})
	add_ginac_test(${tst})
endforeach()

foreach(tmr ${ginac_timings})
	add_ginac_timing(${tmr})
endforeach()

