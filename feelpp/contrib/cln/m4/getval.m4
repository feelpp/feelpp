dnl CLN_HEADER_GETVAL(NAME,FILE)
dnl Expand at autoconf time to the value of a "#define NAME" from the given
dnl FILE. The regexps here aren't very robust, but are enough for us.
dnl /dev/null as a parameter prevents a hang if $2 is accidentally omitted.
dnl (shamelessly ripped from GMP, and changed prefix to CL_).

define(CL_HEADER_GETVAL,
[patsubst(patsubst(
esyscmd([grep "^#define $1 " $2 /dev/null 2>/dev/null]),
[^.*$1[ 	]+],[]),
[[
 	]*$],[])])


dnl CL_GET_VERSION
dnl The CLN version number, extracted from #defines at autoconf time.

AC_DEFUN([CL_GET_VERSION],
[CL_HEADER_GETVAL(CL_VERSION_$1,[include/cln/version.h])])

dnl CL_GET_LTVERSION
dnl The CLN library version number, extracted from #defines at autoconf time.

AC_DEFUN([CL_GET_LTVERSION],
[CL_HEADER_GETVAL(CL_LT_$1,[include/cln/version.h])])
