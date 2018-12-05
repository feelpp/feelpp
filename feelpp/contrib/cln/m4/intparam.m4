# intparam.m4 serial 1
dnl Copyright (C) 2005 Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

dnl From Bruno Haible.

AC_DEFUN([CL_INTPARAM_CROSS],
[
  AC_REQUIRE([CL_LONGLONG])
  AC_REQUIRE([CL_LONGDOUBLE])
  AC_REQUIRE([AC_C_BIGENDIAN])
  cl_machine_file_h=$1
  {
    CL_INTPARAM_BITSIZE([signed char], [char_bitsize])
    CL_INTPARAM_BITSIZE([short], [short_bitsize])
    CL_INTPARAM_BITSIZE([int], [int_bitsize])
    CL_INTPARAM_BITSIZE([long], [long_bitsize])
    if test $cl_cv_c_longlong = yes; then
      CL_INTPARAM_BITSIZE([long long], [longlong_bitsize])
    fi
    CL_INTPARAM_BITSIZE([unsigned char], [uchar_bitsize])
    CL_INTPARAM_BITSIZE([unsigned short], [ushort_bitsize])
    CL_INTPARAM_BITSIZE([unsigned int], [uint_bitsize])
    CL_INTPARAM_BITSIZE([unsigned long], [ulong_bitsize])
    if test $cl_cv_c_longlong = yes; then
      CL_INTPARAM_BITSIZE([unsigned long long], [ulonglong_bitsize])
    fi
    if test -n "$char_bitsize"; then
      echo "/* Integers of type char have $char_bitsize bits. */"
      echo "#define char_bitsize $char_bitsize"
      echo
    else
      echo "#error \"Integers of type char have no binary representation!!\""
    fi
    if test -n "$short_bitsize"; then
      echo "/* Integers of type short have $short_bitsize bits. */"
      echo "#define short_bitsize $short_bitsize"
      echo
    else
      echo "#error \"Integers of type short have no binary representation!!\""
    fi
    if test -n "$int_bitsize"; then
      echo "/* Integers of type int have $int_bitsize bits. */"
      echo "#define int_bitsize $int_bitsize"
      echo
    else
      echo "#error \"Integers of type int have no binary representation!!\""
    fi
    if test -n "$long_bitsize"; then
      echo "/* Integers of type long have $long_bitsize bits. */"
      echo "#define long_bitsize $long_bitsize"
      echo
    else
      echo "#error \"Integers of type long have no binary representation!!\""
    fi
    if test $cl_cv_c_longlong = yes; then
      if test -n "$longlong_bitsize"; then
        echo "/* Integers of type long long have $longlong_bitsize bits. */"
        echo "#define long_long_bitsize $longlong_bitsize"
        echo
      else
        echo "#error \"Integers of type long long have no binary representation!!\""
      fi
    fi
    if test -n "$uchar_bitsize"; then
      echo "/* Integers of type unsigned char have $uchar_bitsize bits. */"
      echo
    else
      echo "#error \"Integers of type unsigned char have no binary representation!!\""
    fi
    if test -n "$ushort_bitsize"; then
      echo "/* Integers of type unsigned short have $ushort_bitsize bits. */"
      echo
    else
      echo "#error \"Integers of type unsigned short have no binary representation!!\""
    fi
    if test -n "$uint_bitsize"; then
      echo "/* Integers of type unsigned int have $uint_bitsize bits. */"
      echo
    else
      echo "#error \"Integers of type unsigned int have no binary representation!!\""
    fi
    if test -n "$ulong_bitsize"; then
      echo "/* Integers of type unsigned long have $ulong_bitsize bits. */"
      echo
    else
      echo "#error \"Integers of type unsigned long have no binary representation!!\""
    fi
    if test $cl_cv_c_longlong = yes; then
      if test -n "$ulonglong_bitsize"; then
        echo "/* Integers of type unsigned long long have $ulonglong_bitsize bits. */"
        echo
      else
        echo "#error \"Integers of type unsigned long long have no binary representation!!\""
      fi
    fi
    if test "$char_bitsize" != "$uchar_bitsize"; then
      echo "#error \"Integer types char and unsigned char have different sizes!!\""
    fi
    if test "$short_bitsize" != "$ushort_bitsize"; then
      echo "#error \"Integer types short and unsigned short have different sizes!!\""
    fi
    if test "$int_bitsize" != "$uint_bitsize"; then
      echo "#error \"Integer types int and unsigned int have different sizes!!\""
    fi
    if test "$long_bitsize" != "$ulong_bitsize"; then
      echo "#error \"Integer types long and unsigned long have different sizes!!\""
    fi
    if test $cl_cv_c_longlong = yes; then
      if test "$longlong_bitsize" != "$ulonglong_bitsize"; then
        echo "#error \"Integer types long long and unsigned long long have different sizes!!\""
      fi
    fi
    AC_TRY_COMPILE([], [typedef int verify[2*(sizeof(char*)<=sizeof (long))-1];],
      [], [echo "#error \"Type char * does not fit into a long!!\""])
    _AC_COMPUTE_INT([sizeof (char *)], [pointer_size])
    pointer_bitsize=`expr $pointer_size '*' $char_bitsize`
    echo "/* Pointers of type char * have $pointer_bitsize bits. */"
    echo "#define pointer_bitsize $pointer_bitsize"
    echo
    CL_INTPARAM_SIZEOF([char], [sizeof_char])
    CL_INTPARAM_ALIGNOF([char], [alignment_char])
    echo "/* Type char has sizeof = $sizeof_char and alignment = $alignment_char. */"
    echo "#define sizeof_char $sizeof_char"
    echo "#define alignment_char $alignment_char"
    echo
    CL_INTPARAM_SIZEOF([unsigned char], [sizeof_uchar])
    CL_INTPARAM_ALIGNOF([unsigned char], [alignment_uchar])
    echo "/* Type unsigned char has sizeof = $sizeof_uchar and alignment = $alignment_uchar. */"
    echo
    CL_INTPARAM_SIZEOF([short], [sizeof_short])
    CL_INTPARAM_ALIGNOF([short], [alignment_short])
    echo "/* Type short has sizeof = $sizeof_short and alignment = $alignment_short. */"
    echo "#define sizeof_short $sizeof_short"
    echo "#define alignment_short $alignment_short"
    echo
    CL_INTPARAM_SIZEOF([unsigned short], [sizeof_ushort])
    CL_INTPARAM_ALIGNOF([unsigned short], [alignment_ushort])
    echo "/* Type unsigned short has sizeof = $sizeof_ushort and alignment = $alignment_ushort. */"
    echo
    CL_INTPARAM_SIZEOF([int], [sizeof_int])
    CL_INTPARAM_ALIGNOF([int], [alignment_int])
    echo "/* Type int has sizeof = $sizeof_int and alignment = $alignment_int. */"
    echo "#define sizeof_int $sizeof_int"
    echo "#define alignment_int $alignment_int"
    echo
    CL_INTPARAM_SIZEOF([unsigned int], [sizeof_uint])
    CL_INTPARAM_ALIGNOF([unsigned int], [alignment_uint])
    echo "/* Type unsigned int has sizeof = $sizeof_uint and alignment = $alignment_uint. */"
    echo
    CL_INTPARAM_SIZEOF([long], [sizeof_long])
    CL_INTPARAM_ALIGNOF([long], [alignment_long])
    echo "/* Type long has sizeof = $sizeof_long and alignment = $alignment_long. */"
    echo "#define sizeof_long $sizeof_long"
    echo "#define alignment_long $alignment_long"
    echo
    CL_INTPARAM_SIZEOF([unsigned long], [sizeof_ulong])
    CL_INTPARAM_ALIGNOF([unsigned long], [alignment_ulong])
    echo "/* Type unsigned long has sizeof = $sizeof_ulong and alignment = $alignment_ulong. */"
    echo
    if test $cl_cv_c_longlong = yes; then
      CL_INTPARAM_SIZEOF([long long], [sizeof_longlong])
      CL_INTPARAM_ALIGNOF([long long], [alignment_longlong])
      echo "/* Type long long has sizeof = $sizeof_longlong and alignment = $alignment_longlong. */"
      echo "#define sizeof_long_long $sizeof_longlong"
      echo "#define alignment_long_long $alignment_longlong"
      echo
      CL_INTPARAM_SIZEOF([unsigned long long], [sizeof_ulonglong])
      CL_INTPARAM_ALIGNOF([unsigned long long], [alignment_ulonglong])
      echo "/* Type unsigned long long has sizeof = $sizeof_ulonglong and alignment = $alignment_ulonglong. */"
      echo
    fi
    CL_INTPARAM_SIZEOF([float], [sizeof_float])
    CL_INTPARAM_ALIGNOF([float], [alignment_float])
    echo "/* Type float has sizeof = $sizeof_float and alignment = $alignment_float. */"
    echo "#define sizeof_float $sizeof_float"
    echo "#define alignment_float $alignment_float"
    echo
    CL_INTPARAM_SIZEOF([double], [sizeof_double])
    CL_INTPARAM_ALIGNOF([double], [alignment_double])
    echo "/* Type double has sizeof = $sizeof_double and alignment = $alignment_double. */"
    echo "#define sizeof_double $sizeof_double"
    echo "#define alignment_double $alignment_double"
    echo
    if test $cl_cv_c_longdouble = yes; then
      CL_INTPARAM_SIZEOF([long double], [sizeof_longdouble])
      CL_INTPARAM_ALIGNOF([long double], [alignment_longdouble])
      echo "/* Type long double has sizeof = $sizeof_longdouble and alignment = $alignment_longdouble. */"
      echo "#define sizeof_long_double $sizeof_longdouble"
      echo "#define alignment_long_double $alignment_longdouble"
      echo
    fi
    CL_INTPARAM_SIZEOF([char *], [sizeof_char_ptr])
    CL_INTPARAM_ALIGNOF([char *], [alignment_char_ptr])
    echo "/* Type char * has sizeof = $sizeof_char_ptr and alignment = $alignment_char_ptr. */"
    echo
    CL_INTPARAM_SIZEOF([long *], [sizeof_long_ptr])
    CL_INTPARAM_ALIGNOF([long *], [alignment_long_ptr])
    echo "/* Type long * has sizeof = $sizeof_long_ptr and alignment = $alignment_long_ptr. */"
    echo
    CL_INTPARAM_SIZEOF([void (*)(void)], [sizeof_function_ptr])
    CL_INTPARAM_ALIGNOF([void (*)(void)], [alignment_function_ptr])
    echo "/* Type function * has sizeof = $sizeof_function_ptr and alignment = $alignment_function_ptr. */"
    echo
    case $ac_cv_c_bigendian in
      yes)
        echo "/* Type unsigned short is stored BIG-ENDIAN in memory (i.e. like mc68000 or sparc). */"
        echo "#define short_big_endian"
        echo "/* Type unsigned int is stored BIG-ENDIAN in memory (i.e. like mc68000 or sparc). */"
        echo "#define int_big_endian"
        echo "/* Type unsigned long is stored BIG-ENDIAN in memory (i.e. like mc68000 or sparc). */"
        echo "#define long_big_endian"
        if test $cl_cv_c_longlong = yes; then
          echo "/* Type unsigned long long is stored BIG-ENDIAN in memory (i.e. like mc68000 or sparc). */"
          echo "#define long_long_big_endian"
        fi
        ;;
      no)
        echo "/* Type unsigned short is stored LITTLE-ENDIAN in memory (i.e. like Z80 or VAX). */"
        echo "#define short_little_endian"
        echo "/* Type unsigned int is stored LITTLE-ENDIAN in memory (i.e. like Z80 or VAX). */"
        echo "#define int_little_endian"
        echo "/* Type unsigned long is stored LITTLE-ENDIAN in memory (i.e. like Z80 or VAX). */"
        echo "#define long_little_endian"
        if test $cl_cv_c_longlong = yes; then
          echo "/* Type unsigned long long is stored LITTLE-ENDIAN in memory (i.e. like Z80 or VAX). */"
          echo "#define long_long_little_endian"
        fi
        ;;
      *)
        echo "#error \"Type short is stored in memory in an obscure manner!!\""
        echo "#error \"Type int is stored in memory in an obscure manner!!\""
        echo "#error \"Type long is stored in memory in an obscure manner!!\""
        if test $cl_cv_c_longlong = yes; then
          echo "#error \"Type long long is stored in memory in an obscure manner!!\""
        fi
        ;;
    esac
    echo
    case $host_cpu in
      hppa)
        echo "/* Stack grows up. */"
        echo "#define stack_grows_up"
        ;;
      *)
        echo "/* Stack grows down. */"
        echo "#define stack_grows_down"
        ;;
    esac
  } > "$cl_machine_file_h"
])

dnl CL_INTPARAM_BITSIZE(type, variable)
dnl puts into variable the determined bitsize of the type.
AC_DEFUN([CL_INTPARAM_BITSIZE],
[
  n=1; x="($1)2"
  while true; do
    AC_TRY_COMPILE([], [typedef int verify[2*(($1)($x) == 0) - 1];],
      [$2=$n; break;],
      [if test $n = 1000; then $2=; break; fi;])
    n=`expr $n + 1`; x="$x * ($1)2"
  done
])

dnl CL_INTPARAM_SIZEOF(type, variable)
dnl puts into variable the determined size of the type.
AC_DEFUN([CL_INTPARAM_SIZEOF],
[
  _AC_COMPUTE_INT([sizeof($1)], [$2])
])

dnl CL_INTPARAM_ALIGNOF(type, variable)
dnl puts into variable the determined alignment of the type.
AC_DEFUN([CL_INTPARAM_ALIGNOF],
[
  dnl Simplify the guessing by assuming that the alignment is a power of 2.
  n=1
  while true; do
    AC_TRY_COMPILE([
#ifdef __cplusplus
# ifdef __GNUC__
#  define alignof(type)  __alignof__ (type)
# else
   template <class type> struct alignof_helper { char slot1; type slot2; };
#  define alignof(type)  offsetof (alignof_helper<type>, slot2)
# endif
#else
# define alignof(type)  offsetof (struct { char slot1; type slot2; }, slot2)
#endif
], [typedef int verify[2*(alignof($1) == $n) - 1];],
      [$2=$n; break;]
      [if test $n = 0; then $2=; break; fi])
    n=`expr $n '*' 2`
  done
])
