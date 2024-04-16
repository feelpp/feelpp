noinst_PROGRAMS = contfrac e fibonacci legendre lucaslehmer nextprime perfnum
bin_PROGRAMS = pi

dist_man_MANS = pi.1


contfrac_SOURCES = contfrac.cc
contfrac_LDADD = ../src/libcln.la

e_SOURCES = e.cc
e_LDADD = ../src/libcln.la

fibonacci_SOURCES = fibonacci.cc
fibonacci_LDADD = ../src/libcln.la

legendre_SOURCES = legendre.cc
legendre_LDADD = ../src/libcln.la

lucaslehmer_SOURCES = lucaslehmer.cc
lucaslehmer_LDADD = ../src/libcln.la

nextprime_SOURCES = nextprime.cc
nextprime_LDADD = ../src/libcln.la

perfnum_SOURCES = perfnum.cc
perfnum_LDADD = ../src/libcln.la

pi_SOURCES = pi.cc
pi_LDADD = ../src/libcln.la


AM_CPPFLAGS = -I../include -I$(top_srcdir)/include
DEFAULT_INCLUDES = -I.@am__isrc@


# Remove .deps directories created by 'configure'.
# I would say that this left-over is an Automake 1.16.1 bug.
distclean-local:
	rm -rf $(DEPDIR)
