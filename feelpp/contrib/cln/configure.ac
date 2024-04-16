dnl run autoreconf -iv to produce the configure script.

m4_define([cl_version_major], CL_GET_VERSION([MAJOR]))
m4_define([cl_version_minor], CL_GET_VERSION([MINOR]))
m4_define([cl_version_patchlevel], CL_GET_VERSION([PATCHLEVEL]))
m4_define([cl_version], [cl_version_major.cl_version_minor.cl_version_patchlevel])

dnl Libtool's library version information for CLN.
m4_define([cl_lt_current], CL_GET_LTVERSION([CURRENT]))
m4_define([cl_lt_revision], CL_GET_LTVERSION([REVISION]))
m4_define([cl_lt_age], CL_GET_LTVERSION([AGE]))
m4_define([cl_lt_version], [cl_lt_current:cl_lt_revision:cl_lt_age])

AC_INIT([cln],[cl_version])
AC_PREREQ([2.69])
AC_CONFIG_SRCDIR(src/integer/gcd/cl_I_gcd.cc)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(autoconf/cl_config.h include/cln/config.h include/cln/host_cpu.h src/base/cl_base_config.h src/base/cl_gmpconfig.h src/timing/cl_t_config.h)
AM_INIT_AUTOMAKE([foreign 1.10 dist-bzip2 subdir-objects])
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
CL_CANONICAL_HOST
                      dnl sets variables host, host_cpu, host_vendor, host_os
CL_CANONICAL_HOST_CPU
                      dnl DEFS __${host_cpu}__
case $host_os in
  dnl Shared libraries are not supported on some platforms [yet]
  mingw* | windows*)
    AC_DISABLE_SHARED
    ;;
  *)
    AC_ENABLE_SHARED
    ;;
esac
LT_INIT
                      dnl sets variable LIBTOOL

LT_VERSION_INFO=cl_lt_version
AC_SUBST(LT_VERSION_INFO)
CL_VERSION=cl_version
AC_SUBST(CL_VERSION)

dnl           checks for compiler characteristics
dnl
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX([11], [ext], [optional])
dnl
AC_CHECK_ALIGNOF([void*])
AC_DEFINE_UNQUOTED([cl_word_alignment], [$ac_cv_alignof_voidp], [Normal alignment of void*])
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
CL_RUSAGE
                      dnl DEFS HAVE_SYS_RESOURCE_H, HAVE_GETRUSAGE, HAVE_SYS_TIMES_H
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
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile examples/Makefile doc/Makefile benchmarks/Makefile cln.spec cln.pc])
AC_OUTPUT
