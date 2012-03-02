
MOSTLYCLEANFILES =
SUPERLUPROGRAMS =

noinst_PROGRAMS = $(SUPERLUPROGRAMS)

AM_CPPFLAGS = 					\
	@AM_CPPFLAGS@ 				\
	$(DUNE_CPPFLAGS) 			\
	$(DUNEMPIFLAGS)

AM_LDFLAGS =					\
	@AM_LDFLAGS@				\
	$(DUNE_LDFLAGS)				\
	$(DUNEMPILDFLAGS)

LDADD =						\
	$(DUNE_LIBS)				\
	$(DUNEMPILIBS)				\
	@LDADD@

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# TRUNK TARGETS FOLLOW

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

if SUPERLU

MOSTLYCLEANFILES += cgstokes_*.vtu cgstokes_*_timings.txt

SUPERLUPROGRAMS += cgstokes
cgstokes_SOURCES = cgstokes.cc
cgstokes_CPPFLAGS =				\
	$(AM_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)			\
	$(UG_CPPFLAGS)				\
	$(ALUGRID_CPPFLAGS)
cgstokes_LDFLAGS =				\
	$(AM_LDFLAGS)				\
	$(SUPERLU_LDFLAGS)			\
	$(UG_LDFLAGS)				\
	$(ALUGRID_LDFLAGS)
cgstokes_LDADD = \
	$(SUPERLU_LIBS)				\
	$(ALUGRID_LIBS)				\
	$(UG_LIBS)

SUPERLUPROGRAMS += cgstokes-velocityblocked
cgstokes_velocityblocked_SOURCES = cgstokes.cc
cgstokes_velocityblocked_CPPFLAGS = $(cgstokes_CPPFLAGS)	\
	-DCGSTOKES_BLOCK_LEVEL=1
cgstokes_velocityblocked_LDFLAGS = $(cgstokes_LDFLAGS)
cgstokes_velocityblocked_LDADD = $(cgstokes_LDADD)

SUPERLUPROGRAMS += cgstokes-allblocked
cgstokes_allblocked_SOURCES = cgstokes.cc
cgstokes_allblocked_CPPFLAGS = $(cgstokes_CPPFLAGS)	\
	-DCGSTOKES_BLOCK_LEVEL=2
cgstokes_allblocked_LDFLAGS = $(cgstokes_LDFLAGS)
cgstokes_allblocked_LDADD = $(cgstokes_LDADD)

SUPERLUPROGRAMS += cgstokes-macroblocks
cgstokes_macroblocks_SOURCES = cgstokes.cc
cgstokes_macroblocks_CPPFLAGS = $(cgstokes_CPPFLAGS)	\
	-DCGSTOKES_MACROBLOCKS
cgstokes_macroblocks_LDFLAGS = $(cgstokes_LDFLAGS)
cgstokes_macroblocks_LDADD = $(cgstokes_LDADD)

SUPERLUPROGRAMS += cgstokes-velocityblocked-macroblocks
cgstokes_velocityblocked_macroblocks_SOURCES = cgstokes.cc
cgstokes_velocityblocked_macroblocks_CPPFLAGS = $(cgstokes_CPPFLAGS)	\
	-DCGSTOKES_BLOCK_LEVEL=1 -DCGSTOKES_MACROBLOCKS
cgstokes_velocityblocked_macroblocks_LDFLAGS = $(cgstokes_LDFLAGS)
cgstokes_velocityblocked_macroblocks_LDADD = $(cgstokes_LDADD)


endif

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
