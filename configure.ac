AC_INIT(dgs,20140613)

AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR(dgs/dgs_gauss_dp.c)

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_LIBTOOL

AC_PROG_INSTALL

AC_CONFIG_HEADERS(dgs/config.h)

AC_PROG_CC_C99()

if test "$ac_cv_prog_cc_c99" = "no"; then
   AC_MSG_ERROR([C99 support is required but not found.])
fi

# Debugging support
AC_ARG_ENABLE(debug, [  --enable-debug          Enable assert() statements for debugging.])

if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile benchmarks/Makefile])
AC_OUTPUT

