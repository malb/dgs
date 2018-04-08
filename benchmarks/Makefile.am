TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/..
STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)

AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPSRCDIR)/dgs -I$(TOPBUILDDIR)-DNDEBUG
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -no-install -ldgs -lmpfr -lgmp -lm

bin_PROGRAMS = bench_gauss bench_rround

.PHONY: clean dist-clean

bench_gauss_SOURCES=bench_gauss.c bench.c bench.h
bench_rround_SOURCES=bench_rround.c bench.c bench.h

