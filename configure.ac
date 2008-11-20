# -*- autoconf -*-
#
#
AC_INIT([life],
        [0.9.9],
        [life-devel@lists.ljkforge.imag.fr])

LIFE_VERSION_MAJOR=0
LIFE_VERSION_MINOR=9
LIFE_VERSION_MICRO=9
LIFE_VERSION_EXTRA=
LIFE_BUILDID=0
LIFE_REVISION=0
LIFE_VERSION="(((${LIFE_VERSION_MAJOR}) << 16) | ((${LIFE_VERSION_MINOR}) << 9) | (${LIFE_VERSION_MICRO}))"
AC_SUBST(LIFE_VERSION_MAJOR)
AC_SUBST(LIFE_VERSION_MINOR)
AC_SUBST(LIFE_VERSION_MICRO)
AC_SUBST(LIFE_VERSION_EXTRA)
AC_SUBST(LIFE_BUILDID)
AC_SUBST(LIFE_REVISION)
AC_SUBST(LIFE_VERSION)
AC_DEFINE_UNQUOTED(LIFE_VERSION_MAJOR,$LIFE_VERSION_MAJOR,[major version number])
AC_DEFINE_UNQUOTED(LIFE_VERSION_MINOR,$LIFE_VERSION_MINOR,[minor version number])
AC_DEFINE_UNQUOTED(LIFE_VERSION_MICRO,$LIFE_VERSION_MICRO,[micro version number])
AC_DEFINE_UNQUOTED(LIFE_VERSION_EXTRA,$LIFE_VERSION_EXTRA,[extra version information: date, cvs, ...])
AC_DEFINE_UNQUOTED(LIFE_BUILDID,$LIFE_BUILDID,[buildid])
AC_DEFINE_UNQUOTED(LIFE_REVISION,$LIFE_REVISION,[revision])
AC_DEFINE_UNQUOTED(LIFE_VERSION,$LIFE_VERSION,[full blown version])

LIFE_PACKAGE=life


AC_DEFINE_UNQUOTED(LIFE_VERSION_STRING,["${LIFE_VERSION_MAJOR}.${LIFE_VERSION_MINOR}.${LIFE_VERSION_MICRO} (${LIFE_VERSION_EXTRA})"],[Life version string])


PACKAGE=$LIFE_PACKAGE
VERSION=${LIFE_VERSION_MAJOR}.${LIFE_VERSION_MINOR}.${LIFE_VERSION_MICRO}
# VERSION=${LIFE_VERSION_MAJOR}.${LIFE_VERSION_MINOR}.${LIFE_VERSION_MICRO}-${LIFE_VERSION_EXTRA}


AC_SUBST(LIFE_PACKAGE)
AC_SUBST(LIFE_VERSION)

# We force the creation of a cache file so configuration in sub-modules will share
# the same values. This is especially true for lib sub-module that set its relative
# path for all subsequent sub-modules.
test -f "$cache_file" || { cache_file=config.cache; :> $cache_file; }

# We set a special cache variable to inform sub-modules that they are built
# together in a combined environment, ie they are not built a stand-alone packages
AC_CACHE_CHECK([for combined build process], [life_cv_combined_env], [life_cv_combined_env=yes])

# Choose our own configuration directory
AC_CONFIG_AUX_DIR([admin])
AC_SUBST([ac_aux_dir])

# Choose project source directory
AC_CONFIG_SRCDIR([configure.ac])


# Propagate through source files the path of the installation directory
AC_DEFINE_UNQUOTED([INSTALL_PATH], "$prefix", [The path of installation])


dnl What the current system and host
AC_CANONICAL_HOST
AC_CANONICAL_SYSTEM


AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADER(lifeconfig.h)

LIFE_CHECK_COMPILERS

dnl check for openmp
ac_save_cxxflags="$CXXFLAGS"
CXXFLAGS="${CXXFLAGS} -fopenmp"
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING(for openmp support in $CXX)
AC_TRY_RUN(AC_LANG_PROGRAM([
#include <omp.h>
],[
omp_set_num_threads(10);
]),
[
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_OPENMP, 1, [Define if the OPENMP is supported by the compiler.])
],
[
  AC_MSG_RESULT(no)
  CXXFLAGS=$ac_save_cxxflags
])


# Turn off shared libraries during beta-testing, since they
# make the build process take too long.
#LT_INIT([disable-static])
#AC_ENABLE_SHARED
AC_PROG_LIBTOOL


# type checks
AC_CHECK_TYPE(size_t, unsigned)

AC_LANG_CPLUSPLUS

AC_MSG_CHECKING(for size_t == uint)
AC_TRY_RUN(AC_LANG_PROGRAM([
#include <cstdlib>
],[
return !(sizeof(size_t)==sizeof(uint));
]),
[
   AC_MSG_RESULT(yes)
   size_t_same_as_uint=yes
   AC_DEFINE(LIFE_SIZET_SAME_AS_UINT, 1, [sizeof(size_t)==sizeof(uint)])
],
[
  AC_MSG_RESULT(no)
  size_t_same_as_uint=no
])


# pdflatex
AC_CHECK_PROG(PDFLATEX, pdflatex,[pdflatex],[no])

# doxygen
AC_CHECK_PROG(DOXYGEN, doxygen, [doxygen], [no-doxygen] )
AC_CHECK_PROG(LIFE_HAS_DOT, dot, [yes] )

# check for boost
with_boost=yes
AX_BOOST([1.33.1])
#LIBS="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SERIALIZATION_LIB} ${BOOST_DATE_TIME_LIB} $LIBS"
theboostlibs="\$(top_builddir)/boost/libs/mpi/src/libboost_mpi.la ${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} \$(top_builddir)/boost/libs/serialization/src/libboost_serialization.la ${BOOST_DATE_TIME_LIB}"
boostlibs="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SIGNALS_LIB} ${BOOST_DATE_TIME_LIB}"

dnl
dnl checks for some scientific computing libraries
dnl
dnl libm
AC_CHECK_LIB(m,cos,[LIBS="-lm ${LIBS}"],${LIBS})
dnl umfpack
# disable it for not since it is not used yet: need to check for "NO"
# AC_CHECK_UMFPACK(4.3.0) # at least umfpack 4.3
AM_CONDITIONAL(HAVE_UMFPACK, test "x$umfpack_libs" != "x" )

dnl check for arbitrary/double-double/quad-double precision types
AC_CHECK_MPFR
AM_CONDITIONAL(HAVE_MPFR, test "x$enable_mpfr" = "xyes")
AC_CHECK_ARPREC
AC_CHECK_QD

dnl lapack, blas, g2c
AC_CHECK_LAPACK

# Check for superlu library
AC_ARG_WITH(superlu-lib,
  [AS_HELP_STRING([--with-superlu-lib],
     [location of the superlu library])],
  with_superlu_lib=$withval, with_superlu_lib="/usr/lib")

AC_CHECK_HEADERS([superlu/slu_Cnames.h superlu/Cnames.h slu_Cnames.h Cnames.h], [CPPFLAGS="${CPPFLAGS} -DGMM_USES_SUPERLU"])

ac_save_ldflags="${LDFLAGS}"
if test "${with_superlu_lib}" != "/usr/lib"; then
 LDFLAGS="-L${with_superlu_lib} ${LDFLAGS}"
fi
SUPERLU_LIBS=
AC_SUBST(SUPERLU_LIBS)
AC_CHECK_LIB(superlu, superlu_abort_and_exit,
             [LIBS="-lsuperlu ${LIBS}";
              with_superlu=yes;
              AC_DEFINE(HAVE_LIBSUPERLU, 1, [Define if the superlu library can be used.])
              ],
             [with_superlu=no; LDFLAGS=$ac_save_ldflags])


# Check for AMD library
AMD_LIBS=
AC_SUBST(AMD_LIBS)
AC_CHECK_LIB(amd, amd_postorder, [AMD_LIBS="-lamd"; with_amd=yes; ],[with_amd=no])

# Check for CAMD library
CAMD_LIBS=
AC_SUBST(CAMD_LIBS)
AC_CHECK_LIB(camd, camd_postorder, [CAMD_LIBS="-lcamd"; with_camd=yes],[with_camd=no])

# Check for UMFPACK library.

UMFPACK_LIBS=
AC_SUBST(UMFPACK_LIBS)

AC_ARG_WITH(umfpack,
  [AS_HELP_STRING([--without-umfpack],
     [don't use UMFPACK, disable some sparse functionality])],
  with_umfpack=$withval, with_umfpack=yes)

warn_umfpack="UMFPACK not found.  This will result in some lack of functionality for sparse matrices."
if test "$with_umfpack" = yes && test "$with_amd" = yes; then
  with_umfpack=no
  AC_CHECK_HEADERS([suitesparse/umfpack.h ufsparse/umfpack.h umfpack/umfpack.h umfpack.h], [
    AC_CHECK_LIB(umfpack, umfpack_zi_get_determinant, [
      UMFPACK_LIBS="-lumfpack"; with_umfpack=yes], [
        ## Invalidate the cache.
	$as_unset ac_cv_lib_umfpack_umfpack_zi_get_determinant
	AC_CHECK_LIB(umfpack, umfpack_zi_get_determinant, [
          UMFPACK_LIBS="-lumfpack"; with_umfpack=yes], [

	  ## Invalidate the cache.
	  $as_unset ac_cv_lib_umfpack_umfpack_zi_get_determinant
	  AC_CHECK_LIB(umfpack, umfpack_zi_get_determinant, [
            UMFPACK_LIBS="-lumfpack -lcblas"; with_umfpack=yes], [], $AMD_LIBS -lcblas $BLAS_LIBS)], $AMD_LIBS $BLAS_LIBS $FLIBS)], $AMD_LIBS)

    if test "$with_umfpack" = yes; then
      AC_DEFINE(HAVE_UMFPACK, 1, [Define if the UMFPACK library is used.])
#      OLD_LIBS=$LIBS
      LIBS="$LIBS $UMFPACK_LIBS $AMD_LIBS $BLAS_LIBS $FLIBS"
#      OCTAVE_UMFPACK_SEPERATE_SPLIT
#      LIBS=$OLD_LIBS
#      TEXINFO_UMFPACK="@set HAVE_UMFPACK"
      warn_umfpack=
    fi
    break])
fi
if test -n "$warn_umfpack"; then
  AC_MSG_WARN($warn_umfpack)
fi

AM_CONDITIONAL(HAVE_UMFPACK, test ! -n "$warn_umfpack"  )

AC_CHECK_MPI
dnl parmetis
AC_CHECK_METIS
AC_CHECK_PARMETIS


dnl petsc
AC_CHECK_PETSC(2.1.6) #
AM_CONDITIONAL(HAVE_PETSC, test "x$petsc_libs" != "x" )


dnl -------------------------------------------------------------
dnl SLEPc  -- disabled by default
dnl -------------------------------------------------------------
AC_ARG_ENABLE(slepc,
              AC_HELP_STRING([--enable-slepc],
                             [build with SLEPc eigen solver support]),
     enableslepc=$enableval,
     enableslepc=yes)

if test "$enableslepc" != no ; then
   LIFE_CHECK_SLEPC
fi
AM_CONDITIONAL(HAVE_SLEPC, test "x$enableslepc" != "xno" )

dnl -------------------------------------------------------------


LIFE_CHECK_TRILINOS


dnl vtk
 AC_CHECK_VTK
#AM_OPTIONS_VTK
#AM_PATH_VTK
#CPPFLAGS=" ${CPPFLAGS} -I/usr/include/vtk-5.0"
#LIBS="-lvtkCommon -lvtkDICOMParser -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkpng -lvtkRendering -lvtktiff -lvtkzlib ${LIBS}"

dnl dnl TAU
dnl ac_save_cppflags="${CPPFLAGS}"
dnl CPPFLAGS="${CPPFLAGS} -I/usr/lib/tau/include"
dnl AC_CHECK_HEADERS([Profile/Profiler.h],[
dnl CPPFLAGS="${CPPFLAGS} \
dnl           -DPROFILING_ON -DTAU_STDCXXLIB -DTAU_GNU \
dnl           -DTAU_DOT_H_LESS_HEADERS -DPTHREADS -DTAU_LARGEFILE\
dnl           -D_LARGEFILE64_SOURCE"],
dnl [CPPFLAGS=$ac_save_cppflags])
dnl ac_save_ldflags="${LDFLAGS}"
dnl LDFLAGS="${LDFLAGS} -L/usr/lib/tau/`uname -m`/lib"
dnl AC_CHECK_LIB(tau-pthread,tau_profile_init,
dnl  [
dnl   AC_DEFINE(HAVE_TAU,1,[defined if tau is installed])
dnl   LIBS="-ltau-pthread ${LIBS}"
dnl  ],
dnl  [LDFLAGS="$ac_save_ldflags"])


dnl dnl grace
dnl AC_CHECK_HEADERS([grace_np.h])
dnl AC_CHECK_LIB([grace_np],GraceClose,[grace_libs="-lgrace_np"],[grace_libs=""])
dnl AC_SUBST(grace_libs)

dnl gmsh
AC_CHECK_GMSH



CPPFLAGS="-I\$(top_srcdir) -I\$(top_srcdir)/boost ${CPPFLAGS} -I\$(top_builddir)"
LIBS="${boostlibs} ${LIBS}"
life_libs="$ac_pwd/life/liblife.la"

AC_MSG_CHECKING([life libs: ])
AC_MSG_RESULT([$life_libs])


AC_SUBST(life_libs)
LIFE_LIBS="-llifefilters -llifediscr -llifepoly -llifemesh  -llifealg -llifecore"

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


#
# setting prefix
#
prefix=${prefix:-"/usr"}
if test "x${prefix}" = "xNONE"; then
    prefix=/usr
fi


##############################################
## Enhance maintainer mode and SUBST variables
## (must be placed after any compilation tests since our overprotective flags
##  let some tests fail)

AM_MAINTAINER_MODE
if test x$USE_MAINTAINER_MODE = xyes
then
   # Maintainers have no choice ! I'm a BOFH, and I plainly assume. [Mt]
   enable_compile_warnings=yes
fi

dnl testdir
AC_CONFIG_TESTDIR(.)

dnl pkgconfig
dnl AX_CREATE_PKGCONFIG_INFO(life.pc,,[${LIFE_LIBS} $LIBS],[Finite Element Library and more],[$CPPFLAGS $CFLAGS])

AC_CONFIG_SUBDIRS([contrib/gmm])
CPPFLAGS="${CPPFLAGS} -I\$(top_srcdir)/contrib/gmm -I\$(top_srcdir)/contrib/gmm/include -I/usr/include/superlu  "

AC_ARG_ENABLE(instantiation-mode,
             AC_HELP_STRING([--enable-instantiation-mode],
                            [build with instantiation support]),
      instantiation_mode=$enableval,
      instantiation_mode=no)

if test x$instantiation_mode = xyes
then
 AC_MSG_NOTICE([build Life with instantiation classes support])
 AC_DEFINE(LIFE_INSTANTIATION_MODE, [1], [Define if support for c++ class instantiation is enabled)])
else
 AC_MSG_NOTICE([build Life without instantiation classes support])
fi

AM_CONDITIONAL([ENABLE_INSTANTIATION_MODE], [test $instantiation_mode = yes])

AC_DEFUN([LIFE_WITH_MESH_MAX_ORDER],
[
 AC_ARG_WITH(meshmaxorder,
             AC_HELP_STRING([--with-meshmaxorder@<:@=INT@:>@],
                            [build with maximum mesh (order @<:@]m4_default([$1], 1)[@:>@ mesh in 1D,2D,3D) support]),
      mesh_max_order=$withval,
      mesh_max_order=1)

 AC_MSG_NOTICE([enabling mesh geometries up to order $mesh_max_order])
 AC_DEFINE(LIFE_MESH_MAX_ORDER, m4_default([$mesh_max_order], 1), [Define if minimal mesh support is required (order 1 mesh in 1D,2D,3D)])

])
LIFE_WITH_MESH_MAX_ORDER

# Activate sub-modules
LIFE_SUBMODULE([doc])
LIFE_SUBMODULE([examples])
LIFE_SUBMODULE([benchmarks])
LIFE_SUBMODULE([applications])


LIFE_CREATE_SUBDIRSLIST
AC_CONFIG_FILES([ Makefile ])


#
# conditional directories
#
AC_ARG_ENABLE([lifesuite],
     [  --enable-lifesuite    Turn on lifesuite creation],
     [case "${enableval}" in
       yes) enable_lifesuite=true ;;
       no)  enable_lifesuite=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-lifesuite]) ;;
     esac],[enable_lifesuite=true])
AM_CONDITIONAL([ENABLE_LIFESUITE], [test $enable_lifesuite = true])

# We force the values of Life install path, headers and library flags
# for subsequent invocations of LIFE_CHECK_LIFE in sub-modules
life_toplevel=
life_install_path=$prefix
life_cppflags="-I\$(top_srcdir)/../ -I\$(top_srcdir)/../boost -I\$(top_builddir)/../ -I\$(top_srcdir)/../../ -I\$(top_srcdir)/../../boost -I\$(top_builddir)/../../  ${CPPFLAGS} "
life_ldflags="${LDFLAGS}"
life_libs="${life_libs} ${LIBS}"

LIFE_SET_LIFE_CACHE_VALUES

#
AC_CONFIG_FILES([ atlocal ])
AC_CONFIG_FILES([ admin/Makefile ])
# boost
AC_CONFIG_FILES([ boost/Makefile ])
AC_CONFIG_FILES([ boost/libs/Makefile ])
AC_CONFIG_FILES([ boost/libs/detail/Makefile ])
AC_CONFIG_FILES([ boost/libs/log/Makefile ])
AC_CONFIG_FILES([ boost/libs/log/src/Makefile ])
AC_CONFIG_FILES([ boost/libs/mpi/Makefile ])
AC_CONFIG_FILES([ boost/libs/mpi/src/Makefile ])
AC_CONFIG_FILES([ boost/libs/serialization/Makefile ])
AC_CONFIG_FILES([ boost/libs/serialization/src/Makefile ])

# contrib
AC_CONFIG_FILES([ contrib/Makefile ])


# life library
AC_CONFIG_FILES([ life/Makefile ])
AC_CONFIG_FILES([ life/lifecore/Makefile ])
AC_CONFIG_FILES([ life/lifealg/Makefile ])
AC_CONFIG_FILES([ life/lifemesh/Makefile ])
AC_CONFIG_FILES([ life/lifepoly/Makefile ])
AC_CONFIG_FILES([ life/lifediscr/Makefile ])
AC_CONFIG_FILES([ life/lifevf/Makefile ])
AC_CONFIG_FILES([ life/lifefilters/Makefile ])
AC_CONFIG_FILES([ life/lifematerial/Makefile ])

# testsuite
AC_CONFIG_FILES([ testsuite/Makefile ])
AC_CONFIG_FILES([ testsuite/lifecore/Makefile ])
AC_CONFIG_FILES([ testsuite/lifealg/Makefile ])
AC_CONFIG_FILES([ testsuite/lifemesh/Makefile ])
AC_CONFIG_FILES([ testsuite/lifepoly/Makefile ])
AC_CONFIG_FILES([ testsuite/lifediscr/Makefile ])
AC_CONFIG_FILES([ testsuite/lifefilters/Makefile ])

AC_DEFINE_UNQUOTED(LIFE_PREFIX, "${prefix:-\"/usr\"}",[define prefix dir for life installation])
AC_OUTPUT


test "x$ENABLE_MODULE_DOC" = x1                  && sum_doc=OK             || sum_doc="no"
test "x$ENABLE_MODULE_EXAMPLES" = x1             && sum_examples=OK        || sum_examples="no"
test "x$ENABLE_MODULE_BENCHMARKS" = x1           && sum_benchmarks=OK      || sum_benchmarks="no"
test "x$ENABLE_MODULE_APPLICATIONS" = x1         && sum_applications=OK    || sum_applications="no"

cat <<EOT

Summary:
========

 Modules              Activation
 -------------------------------
 Doc                  ${sum_doc}
 Examples             ${sum_examples}
 Benchmarks           ${sum_benchmarks}
 Applications         ${sum_applications}

Here are the submodule specific summaries:

EOT

for module in doc examples benchmarks applications
do
  test -f $module/summary && cat $module/summary
done

cat << EOT

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
