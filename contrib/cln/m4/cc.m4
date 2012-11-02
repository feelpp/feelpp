dnl Checks whether the stack can be marked nonexecutable by passing an option
dnl to the C-compiler when acting on .s files. Appends that option to ASFLAGS.
dnl This macro is adapted from one found in GLIBC-2.3.5.
AC_DEFUN([CL_AS_NOEXECSTACK],[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AM_PROG_AS])
AC_CACHE_CHECK([whether --noexecstack is desirable for .s files], cl_cv_as_noexecstack, [dnl
  cat > conftest.c <<EOF
void foo() {}
EOF
  if AC_TRY_COMMAND([${CCAS} $CFLAGS $CPPFLAGS
                     -S -o conftest.s conftest.c >/dev/null]) \
     && grep -q .note.GNU-stack conftest.s \
     && AC_TRY_COMMAND([${CCAS} $CFLAGS $CPPFLAGS -Wa,--noexecstack
                       -c -o conftest.o conftest.s >/dev/null])
  then
    cl_cv_as_noexecstack=yes
  else
    cl_cv_as_noexecstack=no
  fi
  rm -f conftest*])
  if test "$cl_cv_as_noexecstack" = yes; then
    CCASFLAGS="$CCASFLAGS -Wa,--noexecstack"
  fi
])


dnl Checks whether the compiler supports __attribute__((flatten)).
AC_DEFUN([CL_ATTRIBUTE_FLATTEN],[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK([whether the compiler supports __attribute__((flatten))], cl_cv_have_attr_flatten, [dnl
  cat > conftest.cc <<EOF
void f() __attribute__((flatten));
EOF
if AC_TRY_COMMAND(${CXX-g++} $CXXFLAGS -c conftest.cc >/dev/null 2>conftest.stderr)
then
  if grep -i "warning" conftest.stderr > /dev/null; then
    cl_cv_have_attr_flatten=no
  else
    cl_cv_have_attr_flatten=yes
  fi
else
  cl_cv_have_attr_flatten=no
fi
rm -f conftest*
])
if test $cl_cv_have_attr_flatten = yes; then
  AC_DEFINE(CL_HAVE_ATTRIBUTE_FLATTEN, ,[Define if compiler supports __attribute__((flatten))])
fi
])
