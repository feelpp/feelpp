prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_INSTALL_PREFIX@
libdir=@CMAKE_INSTALL_FULL_LIBDIR@
includedir=@CMAKE_INSTALL_FULL_INCLUDEDIR@

Name: cln
Description: Class Library for Numbers
Version: @CL_VERSION@
Libs: @CLN_PC_RPATH@ -L${libdir} -lcln
Libs.private: @GMP_LIBDIR_PC@ -lgmp
Cflags: -I${includedir}
