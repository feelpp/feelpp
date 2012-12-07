prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@EXEC_INSTALL_PREFIX@
libdir=@LIB_INSTALL_DIR@
includedir=@INCLUDE_INSTALL_DIR@

Name: GiNaC
Description: C++ library for symbolic calculations
Version: @GINAC_VERSION@
Requires: cln >= 1.2.2
Libs: -L${libdir} -lginac @GINACLIB_RPATH@
Cflags: -I${includedir}
