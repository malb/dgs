AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${DEBUG_CFLAGS}

lib_LTLIBRARIES = libdgs.la

libdgs_la_SOURCES = dgs/dgs_bern.c \
	dgs/dgs_gauss_dp.c	\
	dgs/dgs_gauss_mp.c	

pkgincludesubdir = $(includedir)/dgs

pkgincludesub_HEADERS = dgs/dgs_bern.h \
	dgs/dgs_gauss.h \
	dgs/dgs.h \
	dgs/dgs_misc.h

libdgs_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined

check_PROGRAMS=test_gauss_z
test_gauss_z_SOURCES=tests/test_gauss_z.c
test_gauss_z_LDFLAGS=-ldgs -lgmp -lmpfr -lm
test_gauss_z_CFLAGS=$(AM_CFLAGS)

TESTS = test_gauss_z

# benchmarketing

SUBDIRS = . benchmarks

clean-local:
	(cd tests; make clean; cd ..)
