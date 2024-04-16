# floatparam.m4 serial 3  -*- Autoconf -*-
dnl Copyright (C) 2005-2008, 2017 Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

dnl From Bruno Haible, Sam Steingold.

AC_DEFUN([CL_FLOATPARAM_CROSS],
[
  cl_machine_file_h=$1
  {
    echo "/* Rounding modes, for use below */"
    echo "#define rounds_to_nearest        0  /* 0.5 ulp */"
    echo "#define rounds_to_zero           1  /* 1 ulp */"
    echo "#define rounds_to_infinity       2  /* 1 ulp */"
    echo "#define rounds_to_minus_infinity 3  /* 1 ulp */"
    echo
    for type in float double "long double"; do
      if test -n "$type"; then
        epsilon_bits=-1; y="($type)1.0"
        while true; do
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  (($type)(($type)1.0 + ($type)($y)) == ($type)1.0)
                  || ($type)(($type)(($type)1.0 + ($type)($y)) - ($type)1.0) != ($type)($y)
                 ) - 1];]])
            ],
            [break;])
          epsilon_bits=`expr $epsilon_bits + 1`; y="$y * ($type)0.5"
        done
        negepsilon_bits=-1; y="($type)-1.0"
        while true; do
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  (($type)(($type)1.0 + ($type)($y)) == ($type)1.0)
                  || ($type)(($type)(($type)1.0 + ($type)($y)) - ($type)1.0) != ($type)($y)
                 ) - 1];]])
            ],
            [break;])
          negepsilon_bits=`expr $negepsilon_bits + 1`; y="$y * ($type)0.5"
        done
        echo "/* Properties of type \`$type': */"
        echo "/* Largest n for which 1+2^(-n) is exactly represented is $epsilon_bits. */"
        echo "/* Largest n for which 1-2^(-n) is exactly represented is $negepsilon_bits. */"
        if test `expr $negepsilon_bits '<=' $epsilon_bits` = 1; then
          echo "#error \"No exponent jump at 1.0 for type $type!\""
        else
          if test `expr $negepsilon_bits '>' $epsilon_bits + 1` = 1; then
            echo "/* Base for type '$type' is 2^"`expr $negepsilon_bits - $epsilon_bits`
          fi
          echo "#define "`echo $type | sed -e 's, ,_,g'`"_mant_bits "`expr $epsilon_bits + 1`
        fi
        x="($type)1.0"
        i=$epsilon_bits
        while test $i != 0; do
          x="$x * ($type)0.5"
          i=`expr $i - 1`
        done
        x="($type)($x)"
        y1="($type)(($type)1.0 + ($type)5.0*$x)"
        y2="($type)(($type)1.0 + ($type)6.0*$x)"
        ys1="($type)(($type)1.0 + ($type)5.4*$x)"
        ys2="($type)(($type)1.0 + ($type)5.6*$x)"
        z1="($type)(($type)-1.0 + ($type)(-5.0)*$x)"
        z2="($type)(($type)-1.0 + ($type)(-6.0)*$x)"
        zs1="($type)(($type)-1.0 + ($type)(-5.4)*$x)"
        zs2="($type)(($type)-1.0 + ($type)(-5.6)*$x)"
        rounds=
        if test -z "$rounds"; then
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  $ys1 == $y1 && $ys2 == $y2 && $zs1 == $z1 && $zs2 == $z2
                 ) - 1];]])
            ],
            [rounds=rounds_to_nearest])
        fi
        if test -z "$rounds"; then
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  $ys1 == $y1 && $ys2 == $y1 && $zs1 == $z1 && $zs2 == $z1
                 ) - 1];]])
            ],
            [rounds=rounds_to_zero])
        fi
        if test -z "$rounds"; then
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  $ys1 == $y2 && $ys2 == $y2 && $zs1 == $z1 && $zs2 == $z1
                 ) - 1];]])
            ],
            [rounds=rounds_to_infinity])
        fi
        if test -z "$rounds"; then
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],
               [[typedef int verify[2*(
                  $ys1 == $y1 && $ys2 == $y1 && $zs1 == $z2 && $zs2 == $z2
                 ) - 1];]])
            ],
            [rounds=rounds_to_minus_infinity])
        fi
        if test -n "$rounds"; then
          echo "#define "`echo $type | sed -e 's, ,_,g'`"_rounds $rounds"
        else
          echo "#error \"Unknown rounding mode for type $type!\""
        fi
        echo
      fi
    done
    dnl Words-in-a-double endianness test. Note that, assuming IEEE 754 format,
    dnl 2.5479915693083957     = { 0x40 0x04 0x62 0x49 0x67 0x65 0x4E 0x64 } ..bIgeNd
    dnl 1.4396527506122064e164 = { 0x62 0x04 0x00 0x00 0x4E 0x65 0x67 0x49 } b...NegI
    dnl 2.5495230282078065     = { 0x40 0x04 0x65 0x6C 0x54 0x54 0x69 0x4C } ..elTTiL
    dnl 1.4139248369879473e214 = { 0x6C 0x65 0x00 0x00 0x4C 0x69 0x54 0x54 } le..LiTT
    double_wordorder_bigendian_p=
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[
         double a[9] = {
           0, 2.5479915693083957, 0, 1.4396527506122064e164,
           0, 2.5495230282078065, 0, 1.4139248369879473e214,
           0 };
         ]],
         [])
      ],
      [if grep LiTTle conftest.$ac_objext >/dev/null; then
         double_wordorder_bigendian_p=0
       else
         if grep bIgeN conftest.$ac_objext >/dev/null; then
           double_wordorder_bigendian_p=1
         fi
       fi
      ])
    if test -n "$double_wordorder_bigendian_p"; then
      echo "#define double_wordorder_bigendian_p $double_wordorder_bigendian_p"
    else
      echo "/* Dazed and confused!  Better not define anything. */"
    fi
    echo
  } > "$cl_machine_file_h"
])
