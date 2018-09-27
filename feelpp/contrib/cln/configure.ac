dnl run autoreconf -iv to produce the configure script.

m4_define([cl_version_major], [1])
m4_define([cl_version_minor], [3])
m4_define([cl_version_patchlevel], [2])
m4_define([cl_version], [cl_version_major.cl_version_minor.cl_version_patchlevel])

dnl Libtool's library version information for CLN.
dnl (Not to be confused with CLN's release version.)
dnl Making new releases:
dnl * increment cl_lt_revision,
dnl * if any interfaces have been added, removed, or changed, then increment
dnl   cl_lt_current and set cl_lt_revision to 0,
dnl * if any interfaces have been added, then increment cl_lt_age,
dnl * if any interfaces have been removed, set cl_lt_age to 0.
dnl (On many systems, $(cl_lt_current):$(cl_lt_revision):$(cl_lt_age) results in
dnl libcln.so.$(cl_lt_current)-$(cl_lt_age).)
m4_define([cl_lt_current], [6])
m4_define([cl_lt_revision], [2])
m4_define([cl_lt_age], [0])
m4_define([cl_lt_version], [cl_lt_current:cl_lt_revision:cl_lt_age])

AC_INIT([cln], cl_version)
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR(src/integer/gcd/cl_I_gcd.cc)
AC_CONFIG_AUX_DIR([autoconf])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(autoconf/cl_config.h include/cln/config.h include/cln/host_cpu.h include/cln/version.h src/base/cl_base_config.h src/base/cl_gmpconfig.h src/timing/cl_t_config.h)
AM_INIT_AUTOMAKE([foreign 1.10 dist-bzip2])
AH_TOP([#ifndef _CL_CONFIG_H])
AH_TOP([#define _CL_CONFIG_H])
AH_TOP([#include "cln/host_cpu.h"])
AH_TOP([/* prevents cln/config.h from being included, so no macro gets redefined */])
AH_TOP([#define _CL_CONFIG_PUBLIC_H])
AH_TOP([/* prevents cln/version.h from being included, so no macro gets redefined */])
AH_TOP([#define _CL_VERSION_H])
AH_BOTTOM([#endif /* _CL_CONFIG_H */])

dnl           checks for programs
dnl
AC_PROG_CC
                      dnl sets variable CC
AC_PROG_CPP
                      dnl sets variable CPP
AC_PROG_CXX
                      dnl sets variable CXX
AC_PROG_CXXCPP
                      dnl sets variable CXXCPP
CL_AS_UNDERSCORE
                      dnl sets variable AS_UNDERSCORE, DEFS ASM_UNDERSCORE
CL_AS_NOEXECSTACK
                      dnl sets variable ASMFLAGS
AC_PROG_INSTALL
                      dnl sets variables INSTALL, INSTALL_DATA, INSTALL_PROGRAM
dnl
dnl           check for host type
dnl
CL_CANONICAL_HOST(autoconf)
                      dnl sets variables host, host_cpu, host_vendor, host_os
CL_CANONICAL_HOST_CPU
                      dnl DEFS __${host_cpu}__
case $host_os in
    *mingw*)
     AC_DISABLE_SHARED
  dnl Shared libraries are not supported on some platforms [yet]
    ;;
    *)
     AC_ENABLE_SHARED
    ;;
esac
AC_PROG_LIBTOOL
                      dnl sets variable LIBTOOL

LT_VERSION_INFO=cl_lt_version
AC_SUBST(LT_VERSION_INFO)

dnl release version
CL_VERSION_MAJOR=cl_version_major
CL_VERSION_MINOR=cl_version_minor
CL_VERSION_PATCHLEVEL=cl_version_patchlevel
dnl release version for cln/config.h, so it can be tested by the preprocessor
AC_DEFINE_UNQUOTED(CL_VERSION_MAJOR, $CL_VERSION_MAJOR, [Major version number of CLN])
AC_DEFINE_UNQUOTED(CL_VERSION_MINOR, $CL_VERSION_MINOR, [Minor version number of CLN])
AC_DEFINE_UNQUOTED(CL_VERSION_PATCHLEVEL, $CL_VERSION_PATCHLEVEL, [Patchlevel version number of CLN])
dnl concatenated release version
CL_VERSION=cl_version
AC_DEFINE_UNQUOTED(CL_VERSION, $CL_VERSION, [CLN release number])
AC_SUBST(CL_VERSION)

dnl           checks for compiler characteristics
dnl
AC_LANG([C++])
dnl
CL_LONGLONG
                      dnl DEFS HAVE_LONGLONG
CL_LONGDOUBLE
                      dnl DEFS HAVE_LONGDOUBLE
dnl
dnl           checks for header files
dnl
AC_CHECK_HEADERS([unistd.h])
                      dnl DEFS HAVE_UNISTD_H if exist unistd.h
dnl
dnl           checks for typedefs
dnl
dnl
dnl           checks for functions and declarations
dnl
CL_ALLOCA
                      dnl set variable ALLOCA, DEFS NO_ALLOCA
CL_GETTIMEOFDAY
                      dnl DEFS HAVE_GETTIMEOFDAY, GETTIMEOFDAY_DOTS, GETTIMEOFDAY_TZP_T
ac_cv_func_ftime=no
                      dnl don't use ftime(3) even if it's available
CL_TIMES_CLOCK
                      dnl DEFS HAVE_TIMES_CLOCK
CL_RUSAGE
                      dnl DEFS HAVE_SYS_RESOURCE_H, HAVE_GETRUSAGE, HAVE_SYS_TIMES_H
CL_PERROR
                      dnl DEFS HAVE_PERROR_DECL
CL_ATTRIBUTE_FLATTEN
                      dnl DEFS CL_HAVE_ATTRIBUTE_FLATTEN
dnl
dnl           checks for compiler characteristics
dnl
AC_C_CHAR_UNSIGNED
                      dnl DEFS __CHAR_UNSIGNED__ if char is unsigned

dnl Create a subdirectory in the ${top_builddir} to hold auto-generated
dnl headers. This subdirectory does not exist when ${top_builddir} != ${top_srcdir}
dnl (a.k.a. "VPATH build").
if test ! -d include/cln; then
	AS_MKDIR_P([include/cln])
fi
CL_MACHINE([integer types and behaviour],${srcdir}/autoconf/intparam.c,[CL_INTPARAM_CROSS],include/cln/intparam.h,cl_cv_file_intparam_h)
                      dnl builds include/cln/intparam.h
CL_MACHINE([floating-point types and behaviour],${srcdir}/autoconf/floatparam.c,[CL_FLOATPARAM_CROSS],include/cln/floatparam.h,cl_cv_file_floatparam_h)
                      dnl builds include/cln/floatparam.h
CL_LIBGMP
                      dnl DEFS CL_USE_GMP, GMP_DEMANDS_UINTD_*

CLNLIB_LIBS='-L${libdir} -lcln'
AC_LIB_LINKFLAGS_FROM_LIBS([CLNLIB_RPATH], [$CLNLIB_LIBS])

dnl
dnl           That's it.
dnl
AC_OUTPUT([Makefile src/Makefile tests/Makefile examples/Makefile doc/Makefile benchmarks/Makefile cln.spec cln.pc])
