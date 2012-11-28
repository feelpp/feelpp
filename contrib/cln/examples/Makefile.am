noinst_PROGRAMS = contfrac e fibonacci legendre lucaslehmer nextprime perfnum
bin_PROGRAMS = pi

dist_man_MANS = pi.1


contfrac_SOURCES = contfrac.cc
contfrac_LDADD = ../src/libfeelpp_cln.la

e_SOURCES = e.cc
e_LDADD = ../src/libfeelpp_cln.la

fibonacci_SOURCES = fibonacci.cc
fibonacci_LDADD = ../src/libfeelpp_cln.la

legendre_SOURCES = legendre.cc
legendre_LDADD = ../src/libfeelpp_cln.la

lucaslehmer_SOURCES = lucaslehmer.cc
lucaslehmer_LDADD = ../src/libfeelpp_cln.la

nextprime_SOURCES = nextprime.cc
nextprime_LDADD = ../src/libfeelpp_cln.la

perfnum_SOURCES = perfnum.cc
perfnum_LDADD = ../src/libfeelpp_cln.la

pi_SOURCES = pi.cc
pi_LDADD = ../src/libfeelpp_cln.la


AM_CPPFLAGS = -I../include -I$(top_srcdir)/include
DEFAULT_INCLUDES = -I.@am__isrc@
