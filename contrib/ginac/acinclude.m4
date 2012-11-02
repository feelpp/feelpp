dnl ===========================================================================
dnl Additional macros used to configure GiNaC.  We don't start our own 
dnl additions' names with AC_ but with GINAC_ in order to steer clear of
dnl future trouble.
dnl ===========================================================================

dnl  GINAC_HEADER_GETVAL(NAME,FILE)
dnl  ----------------------------
dnl  Expand at autoconf time to the value of a "#define NAME" from the given
dnl  FILE. The regexps here aren't very rugged, but are enough for us.
dnl  /dev/null as a parameter prevents a hang if $2 is accidentally omitted.
dnl  (shamelessly ripped from GMP, and changed prefix to GINAC_).

define(GINAC_HEADER_GETVAL,
[patsubst(patsubst(
esyscmd([grep "^#define $1 " $2 /dev/null 2>/dev/null]),
[^.*$1[ 	]+],[]),
[[
 	]*$],[])])
define(GINAC_GET_VERSION,
[GINAC_HEADER_GETVAL(GINACLIB_$1_VERSION,[ginac/version.h])])
define(GINAC_GET_LTVERSION,
[GINAC_HEADER_GETVAL(GINAC_LT_$1,[ginac/version.h])])

dnl Usage: GINAC_STD_CXX_HEADERS
dnl Check for standard C++ headers, bail out if something is missing.
AC_DEFUN([GINAC_STD_CXX_HEADERS], [
AC_CACHE_CHECK([for standard C++ header files], [ginac_cv_std_cxx_headers], [
	ginac_cv_std_cxx_headers="no"
	AC_LANG_PUSH([C++])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
		#include <iosfwd>
		#include <iostream>
		#include <vector>
		#include <list>
		#include <map>
		#include <set>
		#include <string>
		#include <sstream>
		#include <typeinfo>
		#include <stdexcept>
		#include <algorithm>
		#include <limits>
		#include <ctime>
		]])], [ginac_cv_std_cxx_headers="yes"])
	AC_LANG_POP([C++])])
if test "${ginac_cv_std_cxx_headers}" != "yes"; then
	AC_MSG_ERROR([Standard ISO C++ 98 headers are missing])
fi
])

dnl Usage: GINAC_LIBREADLINE
dnl
dnl Check if GNU readline library and headers are avialable.
dnl Defines GINSH_LIBS variable, and HAVE_LIBREADLINE,
dnl HAVE_READLINE_READLINE_H, HAVE_READLINE_HISTORY_H preprocessor macros.
dnl
dnl Note: this macro rejects readline versions <= 4.2 and non-GNU
dnl implementations.
dnl
AC_DEFUN([GINAC_READLINE],[
AC_REQUIRE([GINAC_TERMCAP])
GINSH_LIBS=""
AC_CHECK_HEADERS(readline/readline.h readline/history.h)
if test "x${ac_cv_header_readline_readline_h}" != "xyes" -o "x${ac_cv_header_readline_history_h}" != "xyes"; then
	GINAC_WARNING([ginsh will not compile, because readline headers could not be found.])
else
	AC_CACHE_CHECK([for version of libreadline], [ginac_cv_rl_supported], [
		ginac_cv_rl_supported="no"
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
			#include <stdio.h>
			#include <readline/readline.h>
			#if !defined(RL_VERSION_MAJOR) || !defined(RL_VERSION_MINOR)
			#error "Ancient/unsupported version of readline"
			#endif]])],
			[ginac_cv_rl_supported="yes"])])
	if test "x${ginac_cv_rl_supported}" != "xyes"; then
		GINAC_WARNING([ginsh will not compile, because of an unsupported version of readline (<= 4.2 or non-GNU).])
	else
		save_LIBS="$LIBS"
		LIBS="$LIBTERMCAP $LIBS"
		AC_CHECK_LIB(readline, readline)
		if test "x${ac_cv_lib_readline_readline}" != "xyes"; then
			GINAC_WARNING([ginsh will not compile, because libreadline could not be found.])
		fi
		GINSH_LIBS="$LIBS"
		LIBS="$save_LIBS"
	fi
fi
AC_SUBST(GINSH_LIBS)])
	
dnl Usage: GINAC_TERMCAP
dnl libreadline is based on the termcap functions.
dnl Some systems have tgetent(), tgetnum(), tgetstr(), tgetflag(), tputs(),
dnl tgoto() in libc, some have it in libtermcap, some have it in libncurses.
dnl When both libtermcap and libncurses exist, we prefer the latter, because
dnl libtermcap is being phased out.
AC_DEFUN([GINAC_TERMCAP],
[LIBTERMCAP=
case $host_os in
*mingw32*)
 ;; dnl no termcap libraries are necessary (need hacked libreadline)
*)
AC_CHECK_FUNCS(tgetent)
if test "x$ac_cv_func_tgetent" = "xyes"; then
    :
else
    AC_CHECK_LIB(ncurses, tgetent, LIBTERMCAP="-lncurses")
    if test -z "$LIBTERMCAP"; then
        AC_CHECK_LIB(termcap, tgetent, LIBTERMCAP="-ltermcap")
    fi
fi
;;
esac
AC_SUBST(LIBTERMCAP)
])

dnl Usage: GINAC_ERROR(message)
dnl This macro displays the warning "message" and sets the flag ginac_error
dnl to yes.
AC_DEFUN([GINAC_ERROR],[
ginac_error_txt="$ginac_error_txt
** $1
"
ginac_error=yes])

dnl Usage: GINAC_WARNING(message)
dnl This macro displays the warning "message" and sets the flag ginac_warning
dnl to yes.
AC_DEFUN([GINAC_WARNING],[
ginac_warning_txt="$ginac_warning_txt
== $1
"
ginac_warning=yes])

dnl Usage: GINAC_CHECK_ERRORS
dnl (must be put at end of configure.in, because it exits on error)
dnl This macro displays a warning message if GINAC_ERROR or GINAC_WARNING 
dnl has occured previously.
AC_DEFUN([GINAC_CHECK_ERRORS],[
if test "x${ginac_error}" = "xyes"; then
	echo
    echo "**** The following problems have been detected by configure."
    echo "**** Please check the messages below before running \"make\"."
    echo "**** (see the section 'Common Problems' in the INSTALL file)"
    echo "$ginac_error_txt"
    if test "x${ginac_warning_txt}" != "x"; then
        echo "${ginac_warning_txt}"
    fi
    if test "x$cache_file" != "x/dev/null"; then
        echo "deleting cache ${cache_file}"
        rm -f $cache_file
    fi
    exit 1
else 
    if test "x${ginac_warning}" = "xyes"; then
		echo
        echo "=== The following minor problems have been detected by configure."
        echo "=== Please check the messages below before running \"make\"."
        echo "=== (see the section 'Common Problems' in the INSTALL file)"
        echo "$ginac_warning_txt"
    fi
    echo "Configuration of GiNaC $VERSION done. Now type \"make\"."
fi])

AC_DEFUN([GINAC_HAVE_RUSAGE],
[AC_CACHE_CHECK([whether struct rusage is declared in <sys/resource.h>],
ac_cv_have_rusage,
 [AC_TRY_COMPILE([#include <sys/times.h>
                  #include <sys/resource.h>],
                  [struct rusage resUsage;
                   getrusage(RUSAGE_SELF, &resUsage);
                   return 0;],
                 [ac_cv_have_rusage=yes],
                 [ac_cv_have_rusage=no])
])
CONFIG_RUSAGE="no"
if test "$ac_cv_have_rusage" = yes; then
  CONFIG_RUSAGE="yes"
  AC_DEFINE(HAVE_RUSAGE,,[define if struct rusage declared in <sys/resource.h>])
fi
AC_SUBST(CONFIG_RUSAGE)
])

dnl Usage: GINAC_EXCOMPILER
dnl - Checks if dlopen is available
dnl - Allows user to disable GiNaC::compile_ex (e.g. for security reasons)
dnl Defines HAVE_LIBDL preprocessor macro, sets DL_LIBS and CONFIG_EXCOMPILER
dnl variables.
AC_DEFUN([GINAC_EXCOMPILER], [
CONFIG_EXCOMPILER=yes
DL_LIBS=""

AC_ARG_ENABLE([excompiler], 
	[AS_HELP_STRING([--enable-excompiler], [Enable GiNaC::compile_ex (default: yes)])],
	[if test "$enableval" = "no"; then
		CONFIG_EXCOMPILER="no"
	fi],
	[CONFIG_EXCOMPILER="yes"])

case $host_os in
	*mingw32*)
	CONFIG_EXCOMPILER="notsupported"
	;;
	*)
	;;
esac

if test "$CONFIG_EXCOMPILER" = "yes"; then
	AC_CHECK_HEADER([dlfcn.h], [CONFIG_EXCOMPILER="yes"], [CONFIG_EXCOMPILER="no"])
elif test "$CONFIG_EXCOMPILER" = "no"; then
	AC_MSG_NOTICE([GiNaC::compile_ex disabled at user request.])
else
	AC_MSG_NOTICE([GiNaC::compile_ex is not supported on $host_os.])
fi
	
if test "$CONFIG_EXCOMPILER" = "yes"; then
	dnl Some systems (GNU/Linux, Solaris) have dlopen in -ldl, some
	dnl others (OpenBSD) -- in libc
	found_dlopen_lib="no"
	DL_LIBS="-ldl"
	AC_CHECK_LIB(dl, dlopen, [found_dlopen_lib="yes"])
	if test "$found_dlopen_lib" = "no"; then
		DL_LIBS=""
		AC_CHECK_FUNC(dlopen, [found_dlopen_lib="yes"])
	fi
	if test "$found_dlopen_lib" = "no"; then
		CONFIG_EXCOMPILER="no"
		AC_MSG_WARN([Could not found working dlopen(). GiNaC::compile_ex will be disabled.])
	else
		AC_DEFINE(HAVE_LIBDL, 1, [set to 1 if dlopen() works.])
	fi
fi
AC_SUBST(DL_LIBS)
AC_SUBST(CONFIG_EXCOMPILER)])

