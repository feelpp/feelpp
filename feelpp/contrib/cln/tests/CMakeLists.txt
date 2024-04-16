
set(cln_tests exam tests)

set(exam_SOURCES
	exam.h exam.cc
	exam_I.cc exam_I_plus.cc exam_I_minus.cc exam_I_mul.cc
	exam_I_div.cc exam_I_floor.cc
	exam_RA.cc exam_RA_plus.cc exam_RA_minus.cc exam_RA_mul.cc
	exam_RA_div.cc exam_RA_floor.cc
	exam_SF.cc exam_SF_plus.cc exam_SF_minus.cc exam_SF_mul.cc
	exam_SF_div.cc exam_SF_floor.cc
	exam_FF.cc exam_FF_plus.cc exam_FF_minus.cc exam_FF_mul.cc
	exam_FF_div.cc exam_FF_floor.cc
	exam_DF.cc exam_DF_plus.cc exam_DF_minus.cc exam_DF_mul.cc
	exam_DF_div.cc exam_DF_floor.cc
	exam_LF.cc exam_LF_plus.cc exam_LF_minus.cc exam_LF_mul.cc
	exam_LF_div.cc exam_LF_floor.cc
	exam_I_factorial.cc
	exam_I_gcd.cc exam_I_sqrtp.cc test_MI.h test.h
)

set(tests_SOURCES
	test.h tests.cc test_I.cc test_I.h test_I_abs.cc test_I_compare.cc
	test_I_plus.cc test_I_minus.cc test_I_plus1.cc
	test_I_minus1.cc test_I_mul.cc test_I_div.cc
	test_I_gcd.cc test_I_xgcd.cc test_I_ash.cc
	test_I_evenp.cc test_I_oddp.cc test_I_lognot.cc
	test_I_logand.cc test_I_logandc1.cc test_I_logandc2.cc
	test_I_logior.cc test_I_logorc1.cc test_I_logorc2.cc
	test_I_logxor.cc test_I_lognand.cc test_I_lognor.cc
	test_I_logeqv.cc test_I_boole.cc test_I_logbitp.cc
	test_I_logtest.cc test_I_ldb.cc test_I_ldbtest.cc
	test_I_mkf.cc test_I_dpb.cc test_I_dpf.cc
	test_I_logcount.cc test_I_ilength.cc test_I_ord2.cc
	test_I_power2p.cc test_I_isqrt.cc test_I_sqrtp.cc
	test_I_io.cc test_I_GV.cc
	test_MI.h test_MI.cc test_MI_canonhom.cc test_MI_plus.cc
	test_MI_minus.cc test_MI_mul.cc test_MI_recip.cc
	test_MI_div.cc test_MI_expt.cc
	test_nt.h test_nt.cc test_nt_jacobi.cc)
	
macro(add_cln_test thename)
	if ("${${thename}_SOURCES}" STREQUAL "")
		set(${thename}_sources ${thename}.cc ${${thename}_extra_src})
	endif()
	add_executable(${thename} EXCLUDE_FROM_ALL ${${thename}_SOURCES})
	target_link_libraries(${thename} cln::cln)
	target_include_directories(${thename}
		PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/../src
		${CMAKE_CURRENT_BINARY_DIR}/../src
		${CMAKE_CURRENT_BINARY_DIR}/../autoconf
	)
	add_dependencies(test_suite ${thename})
	add_dependencies(check ${thename})
	add_test(NAME ${thename} COMMAND $<TARGET_FILE:${thename}>)
endmacro()

foreach(tst ${cln_tests})
	add_cln_test(${tst})
endforeach()

