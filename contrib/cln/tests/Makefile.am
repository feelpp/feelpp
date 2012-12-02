
ALL_TESTS = \
exam \
tests

ALL_TIMINGS = \
timemul \
timesquare \
timediv \
timesqrt \
timegcd \
timefact \
timeprint \
timeLFsqrt \
timeRAtoLF \
timeLFRAmul \
timeRALFdiv \
timepi \
timeLFln \
timeLFexp \
timeLFsin \
timeLFcos \
timeLFsinh \
timeLFcosh \
timeLFatan \
timeLFatanh \
timeMIpow2recip \
timeMIpow2div \
timeMImisc5 \
timeUPMImul \
timesqrtmodp \
timeexp1 \
timeeuler \
timecatalan \
timezeta3 \
timerecip2adic \
timediv2adic \
main


TESTS = $(ALL_TESTS)
check_PROGRAMS = $(ALL_TESTS)
EXTRA_PROGRAMS = $(ALL_TIMINGS)

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src \
-I$(top_builddir)/include -I$(top_builddir)/src
DEFAULT_INCLUDES = -I.@am__isrc@ -I$(top_builddir)/autoconf


exam_SOURCES = exam.h exam.cc \
	       exam_I.cc exam_I_plus.cc exam_I_minus.cc exam_I_mul.cc \
	       exam_I_div.cc exam_I_floor.cc \
	       exam_RA.cc exam_RA_plus.cc exam_RA_minus.cc exam_RA_mul.cc \
	       exam_RA_div.cc exam_RA_floor.cc \
	       exam_SF.cc exam_SF_plus.cc exam_SF_minus.cc exam_SF_mul.cc \
	       exam_SF_div.cc exam_SF_floor.cc \
	       exam_FF.cc exam_FF_plus.cc exam_FF_minus.cc exam_FF_mul.cc \
	       exam_FF_div.cc exam_FF_floor.cc \
	       exam_DF.cc exam_DF_plus.cc exam_DF_minus.cc exam_DF_mul.cc \
	       exam_DF_div.cc exam_DF_floor.cc \
	       exam_LF.cc exam_LF_plus.cc exam_LF_minus.cc exam_LF_mul.cc \
	       exam_LF_div.cc exam_LF_floor.cc \
	       exam_I_gcd.cc exam_I_sqrtp.cc test_MI.h test.h

exam_LDADD = ../src/libfeelpp_cln.la

tests_SOURCES = test.h tests.cc test_I.cc test_I.h test_I_abs.cc test_I_compare.cc \
	test_I_plus.cc test_I_minus.cc test_I_plus1.cc \
	test_I_minus1.cc test_I_mul.cc test_I_div.cc \
	test_I_gcd.cc test_I_xgcd.cc test_I_ash.cc \
	test_I_evenp.cc test_I_oddp.cc test_I_lognot.cc \
	test_I_logand.cc test_I_logandc1.cc test_I_logandc2.cc \
	test_I_logior.cc test_I_logorc1.cc test_I_logorc2.cc \
	test_I_logxor.cc test_I_lognand.cc test_I_lognor.cc \
	test_I_logeqv.cc test_I_boole.cc test_I_logbitp.cc \
	test_I_logtest.cc test_I_ldb.cc test_I_ldbtest.cc \
	test_I_mkf.cc test_I_dpb.cc test_I_dpf.cc \
	test_I_logcount.cc test_I_ilength.cc test_I_ord2.cc \
	test_I_power2p.cc test_I_isqrt.cc test_I_sqrtp.cc \
	test_I_io.cc test_I_GV.cc \
	test_MI.h test_MI.cc test_MI_canonhom.cc test_MI_plus.cc \
	test_MI_minus.cc test_MI_mul.cc test_MI_recip.cc \
	test_MI_div.cc test_MI_expt.cc \
	test_nt.h test_nt.cc test_nt_jacobi.cc

tests_LDADD = ../src/libfeelpp_cln.la

timemul_SOURCES = timemul.cc
timemul_LDADD = ../src/libfeelpp_cln.la

timesquare_SOURCES = timesquare.cc
timesquare_LDADD = ../src/libfeelpp_cln.la

timediv_SOURCES = timediv.cc
timediv_LDADD = ../src/libfeelpp_cln.la

timesqrt_SOURCES = timesqrt.cc
timesqrt_LDADD = ../src/libfeelpp_cln.la

timegcd_SOURCES = timegcd.cc
timegcd_LDADD = ../src/libfeelpp_cln.la

timefact_SOURCES = timefact.cc
timefact_LDADD = ../src/libfeelpp_cln.la

timeprint_SOURCES = timeprint.cc
timeprint_LDADD = ../src/libfeelpp_cln.la

timeLFsqrt_SOURCES = timeLFsqrt.cc
timeLFsqrt_LDADD = ../src/libfeelpp_cln.la

timeRAtoLF_SOURCES = timeRAtoLF.cc
timeRAtoLF_LDADD = ../src/libfeelpp_cln.la

timeLFRAmul_SOURCES = timeLFRAmul.cc
timeLFRAmul_LDADD = ../src/libfeelpp_cln.la

timeRALFdiv_SOURCES = timeRALFdiv.cc
timeRALFdiv_LDADD = ../src/libfeelpp_cln.la

timepi_SOURCES = timepi.cc
timepi_LDADD = ../src/libfeelpp_cln.la

timeexp1_SOURCES = timeexp1.cc
timeexp1_LDADD = ../src/libfeelpp_cln.la

timeeuler_SOURCES = timeeuler.cc
timeeuler_LDADD = ../src/libfeelpp_cln.la

timecatalan_SOURCES = timecatalan.cc
timecatalan_LDADD = ../src/libfeelpp_cln.la

timezeta3_SOURCES = timezeta3.cc
timezeta3_LDADD = ../src/libfeelpp_cln.la

timeLFln_SOURCES = timeLFln.cc
timeLFln_LDADD = ../src/libfeelpp_cln.la

timeLFexp_SOURCES = timeLFexp.cc
timeLFexp_LDADD = ../src/libfeelpp_cln.la

timeLFsin_SOURCES = timeLFsin.cc
timeLFsin_LDADD = ../src/libfeelpp_cln.la

timeLFcos_SOURCES = timeLFcos.cc
timeLFcos_LDADD = ../src/libfeelpp_cln.la

timeLFsinh_SOURCES = timeLFsinh.cc
timeLFsinh_LDADD = ../src/libfeelpp_cln.la

timeLFcosh_SOURCES = timeLFcosh.cc
timeLFcosh_LDADD = ../src/libfeelpp_cln.la

timeLFatan_SOURCES = timeLFatan.cc
timeLFatan_LDADD = ../src/libfeelpp_cln.la

timeLFatanh_SOURCES = timeLFatanh.cc
timeLFatanh_LDADD = ../src/libfeelpp_cln.la

timerecip2adic_SOURCES = timerecip2adic.cc
timerecip2adic_LDADD = ../src/libfeelpp_cln.la

timediv2adic_SOURCES = timediv2adic.cc
timediv2adic_LDADD = ../src/libfeelpp_cln.la

timeMIpow2recip_SOURCES = timeMIpow2recip.cc
timeMIpow2recip_LDADD = ../src/libfeelpp_cln.la

timeMIpow2div_SOURCES = timeMIpow2div.cc
timeMIpow2div_LDADD = ../src/libfeelpp_cln.la

timeMImisc5_SOURCES = timeMImisc5.cc
timeMImisc5_LDADD = ../src/libfeelpp_cln.la

timeUPMImul_SOURCES = timeUPMImul.cc
timeUPMImul_LDADD = ../src/libfeelpp_cln.la

timesqrtmodp_SOURCES = timesqrtmodp.cc
timesqrtmodp_LDADD = ../src/libfeelpp_cln.la

main_SOURCES = main.cc
main_LDADD = ../src/libfeelpp_cln.la


