# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/course-examples/example01.cc])
AM_CONFIG_HEADER([config.h])

AC_CHECK_PROGS([CONVERT], [convert], [false])
AM_CONDITIONAL([CONVERT], [test "x$CONVERT" != xfalse])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-localfunctions dune-pdelab
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-localfunctions], [dune-pdelab]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  m4/Makefile
  doc/Makefile
  doc/howto/Makefile
  src/Makefile
  src/utility/Makefile
  src/course-examples/Makefile
  src/course-gridfunctionspace/Makefile
  src/navier-stokes/Makefile
  src/two-phase/Makefile
  src/convection-diffusion/Makefile
  dune-pdelab-howto.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
