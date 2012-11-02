SUBDIRS = src tests examples doc benchmarks
DIST_SUBDIRS = src tests examples doc benchmarks

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = cln.pc

EXTRA_DIST = cln.pc.in cln.spec.in
DISTCLEANFILES = cln.pc

$(pkgconfig_DATA): config.status

# Rule to build tar-bzipped distribution package
$(PACKAGE)-$(VERSION).tar.gz: dist

ACLOCAL_AMFLAGS = -I m4


clnincludedir = $(includedir)/cln

clninclude_HEADERS = \
include/cln/SV_real.h \
include/cln/GV_integer.h \
include/cln/floatformat.h \
include/cln/lfloat.h \
include/cln/null_ring.h \
include/cln/io.h \
include/cln/ring.h \
include/cln/V.h \
include/cln/GV.h \
include/cln/GV_number.h \
include/cln/complex_io.h \
include/cln/condition.h \
include/cln/integer_ring.h \
include/cln/dfloat_io.h \
include/cln/random.h \
include/cln/SV_ringelt.h \
include/cln/ffloat_io.h \
include/cln/float_class.h \
include/cln/cln.h \
include/cln/ffloat_class.h \
include/cln/float.h \
include/cln/string.h \
include/cln/lfloat_io.h \
include/cln/malloc.h \
include/cln/lfloat_class.h \
include/cln/number_io.h \
include/cln/numtheory.h \
include/cln/object.h \
include/cln/proplist.h \
include/cln/univpoly_rational.h \
include/cln/univpoly_real.h \
include/cln/output.h \
include/cln/real_io.h \
include/cln/ffloat.h \
include/cln/sfloat_io.h \
include/cln/timing.h \
include/cln/SV_number.h \
include/cln/complex_ring.h \
include/cln/univpoly_complex.h \
include/cln/version.h \
include/cln/rational.h \
include/cln/rational_class.h \
include/cln/rational_io.h \
include/cln/types.h \
include/cln/univpoly_modint.h \
include/cln/modinteger.h \
include/cln/rational_ring.h \
include/cln/univpoly_integer.h \
include/cln/number.h \
include/cln/GV_complex.h \
include/cln/GV_modinteger.h \
include/cln/GV_real.h \
include/cln/SV_complex.h \
include/cln/SV_integer.h \
include/cln/complex.h \
include/cln/exception.h \
include/cln/univpoly.h \
include/cln/SV_rational.h \
include/cln/complex_class.h \
include/cln/real.h \
include/cln/symbol.h \
include/cln/dfloat_class.h \
include/cln/modules.h \
include/cln/real_ring.h \
include/cln/float_io.h \
include/cln/GV_rational.h \
include/cln/input.h \
include/cln/integer_class.h \
include/cln/integer_io.h \
include/cln/real_class.h \
include/cln/sfloat.h \
include/cln/sfloat_class.h \
include/cln/dfloat.h \
include/cln/SV.h \
include/cln/integer.h

## FIXME: we should NOT expose these macros to users
nodist_clninclude_HEADERS = \
include/cln/config.h \
include/cln/host_cpu.h \
include/cln/intparam.h

DISTCLEANFILES += \
include/cln/config.h \
include/cln/host_cpu.h \
include/cln/intparam.h \
include/cln/floatparam.h

EXTRA_DIST += include/cln/config.h.in \
	     include/cln/host_cpu.h.in \
	     include/cln/version.h.in \
autoconf/floatparam.c \
autoconf/intparam.c

DISTCHECK_CONFIGURE_FLAGS = --disable-static

