BENCHCFLAGS = $(AM_CFLAGS)  -I.. @CFLAGS@ -DNDEBUG
BENCHLIBADD = -lmpfr -lgmp -lm
BENCHLDFLAGS = -Wl,-rpath,../.libs/ ../.libs/libdgs.la

bin_PROGRAMS = bench_gauss

.PHONY: clean dist-clean

bench_gauss_SOURCES=bench_gauss.c bench.c bench.h
bench_gauss_CFLAGS=${BENCHCFLAGS}
bench_gauss_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

