#                                               -*- Autoconf -*-
#
#  configure.ac
#
#  (C) Copyright 2008 Université de Grenoble 1 (Joseph Fourier)
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  Process this file with autoconf to produce a configure script.
#  Set minimal version of Autoconf for configure to work
#
AC_PREREQ([2.59])

# Initialize Autoconf
AC_INIT([Life platform module benchmarks],
        [0.9.11],
        [life-commits@lists.ljkforge.imag.fr],
        [life-benchmarks])

# We check a special cache variable to see if we are built in a combined environment
AC_CACHE_CHECK([for combined build process], [life_cv_combined_env], [life_cv_combined_env=no])

# Choose our own configuration directory
AC_CONFIG_AUX_DIR([admin])
AC_SUBST([ac_aux_dir])

# Initialize Automake
AM_INIT_AUTOMAKE([foreign])

# Chooze project source directory
AC_CONFIG_SRCDIR([configure.ac])

# Propagate through source files the path of the installation directory
AC_DEFINE_UNQUOTED([INSTALL_PATH], "$prefix", [The path of installation])


# Proceed to required program detection
AC_LIBTOOL_DLOPEN # MUST be declared BEFORE AC_PROG_LIBTOOL
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS
AC_PROG_CPP
AC_PROG_YACC
AC_PROG_LEX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Automatically save cache file to improve rerunning of configure script
AC_CACHE_SAVE

LIFE_DEBUG
LIFE_CHECK_COMPILERS
LIFE_CHECK_LIFE

dnl ------------------------------------------------------------
dnl Python
dnl ------
LIFE_CHECK_PYTHON
LIFE_CHECK_PYTHON_MODULE([scipy])
LIFE_CHECK_PYTHON_MODULE([pyx])

# Write out configure arguments
life_configure_args=`eval echo configure $ac_configure_args`
AC_SUBST(life_configure_args)

# Write out configuration date
AC_PATH_PROG([DATE], [date], [:])
life_configure_date=`$DATE -R`
AC_SUBST(life_configure_date)

# Write out project revision if versioning tool is present
AC_PATH_PROG([SVNVERSION], [svnversion], [echo unknown])
#_AC_SRCDIRS([.])
life_revision=`cd "$ac_top_srcdir" && $SVNVERSION`
AC_SUBST(life_revision)

# Automatically save cache file to improve rerunning of configure script
AC_CACHE_SAVE


AC_CONFIG_TESTDIR(.)

#
LIFE_EXTRAMODULE([validation])

# benchmark
AC_CONFIG_FILES([ atlocal ])
AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ admin/Makefile ])
AC_CONFIG_FILES([ convergence/Makefile ])
AC_CONFIG_FILES([ perf/Makefile ])

#AC_CONFIG_FILES([ kovasznay/Makefile ])
#AC_CONFIG_FILES([ stokes/Makefile ])
#AC_CONFIG_FILES([ turek/Makefile ])
#AC_CONFIG_FILES([ structure/Makefile ])

AC_OUTPUT

test "x$ENABLE_EXTRAMODULE_VALIDATION" = xvalidation      && sum_benchmarks_validation=OK    || sum_benchmarks_validation="no"

test "x$WITH_PYTHON" = x1                && sum_python=OK                   || sum_python="no * (optional)"
test "x$WITH_PYTHON_MODULE_SCIPY" = x1   && sum_python_module_scipy=OK      || sum_python_module_scipy="no * (optional)"
test "x$WITH_PYTHON_MODULE_PYX" = x1     && sum_python_module_pyx=OK        || sum_python_module_pyx="no * (optional)"

cat <<EOT > summary

Benchmarks Sub-Module Summary:
==============================

 Extra-Modules        Activation
 -------------------------------
 Validation           ${sum_benchmarks_validation}

(lib)
 Debug level is ${DEBUG_LEVEL}

 Support for          is
 -----------------------------

  Python              ${sum_python}
  Python module scipy ${sum_python_module_scipy}
  Python module pyx   ${sum_python_module_pyx}

EOT
cat summary

cat <<EOT

Congratulations !
Now you can build Life as usual:
	make
	make check (optional but recommended)
	make install
	make installcheck (optional but recommended)

And if you find the output too verbose, try using the following options:
	make -s LIBTOOL_FLAGS=--silent

Have fun !
EOT
