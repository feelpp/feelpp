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
AC_INIT([Life platform module documentation],
        [0.9.10],
        [life-commits@lists.ljkforge.imag.fr],
        [life-doc])

# We check a special cache variable to see if we are built in a combined environment
AC_CACHE_CHECK([for combined build process], [life_cv_combined_env], [life_cv_combined_env=no])

# Choose our own configuration directory
AC_CONFIG_AUX_DIR([admin])
AC_SUBST([ac_aux_dir])

# Initialize Automake
AM_INIT_AUTOMAKE([foreign])

AM_MAINTAINER_MODE

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

LIFE_DEBUG
LIFE_CHECK_COMPILERS
LIFE_CHECK_LIFE

# Automatically save cache file to improve rerunning of configure script
AC_CACHE_SAVE

LIFE_CHECK_PDFLATEX
LIFE_CHECK_MAKEINDEX
LIFE_CHECK_BIBTEX

LIFE_CHECK_DOXYGEN
LIFE_CHECK_DOT
LIFE_CHECK_FIG2DEV

# Automatically save cache file to improve rerunning of configure script
AC_CACHE_SAVE

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

# doc
AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ admin/Makefile ])
AC_CONFIG_FILES([ tutorial/Makefile ])
AC_CONFIG_FILES([ tutorial/version.tex ])
AC_CONFIG_FILES([ tutorial/styles/Makefile ])
AC_CONFIG_FILES([ tutorial/pdfs/Makefile ])
AC_CONFIG_FILES([ tutorial/pngs/Makefile ])
AC_CONFIG_FILES([ api/Makefile api/Doxyfile ])

AC_OUTPUT

test "x$WITH_PDFLATEX" = x1              && sum_pdflatex=OK                 || sum_pdflatex="no * (optional)"
test "x$WITH_MAKEINDEX" = x1             && sum_makeindex=OK                || sum_makeindex="no * (optional)"
test "x$WITH_BIBTEX" = x1                && sum_bibtex=OK                   || sum_bibtex="no * (optional)"
test "x$WITH_DOXYGEN" = x1               && sum_doxygen=OK                  || sum_doxygen="no * (optional)"
test "x$WITH_DOT" = x1                   && sum_dot=OK                      || sum_dot="no * (optional)"
test "x$WITH_FIG2DEV" = x1               && sum_fig2dev=OK                  || sum_fig2dev="no * (optional)"


cat <<EOT > summary

(lib)
 Debug level is ${DEBUG_LEVEL}

 Support for          is
 -----------------------------
  PDFLatex            ${sum_pdflatex}
  MakeIndex           ${sum_makeindex}
  BibTeX              ${sum_bibtex}
  Doxygen             ${sum_doxygen}
  Dot                 ${sum_dot}
  Fig2dev             ${sum_fig2dev}

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

