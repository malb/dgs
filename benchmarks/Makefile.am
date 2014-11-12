BENCHCFLAGS = $(AM_CFLAGS)  -I.. @CFLAGS@ -DNDEBUG
BENCHLIBADD = -lgmp -lmpfr -lm
BENCHLDFLAGS = -Wl,-rpath,../.libs/ ../.libs/libdgs.la

EXTRA_PROGRAMS = bench_gauss

all: bench_gauss

.PHONY: clean dist-clean

clean-local:
	-rm -f $(EXTRA_PROGRAMS)

bench_gauss_SOURCES=bench_gauss.c bench.c
bench_gauss_CFLAGS=${BENCHCFLAGS}
bench_gauss_LDFLAGS=${BENCHLDFLAGS} ${BENCHLIBADD}

