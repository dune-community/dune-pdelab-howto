
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

if SUPERLU

MOSTLYCLEANFILES += cgstokes_*.vtu cgstokes_*_timings.txt

SUPERLUPROGRAMS += cgstokes-old
cgstokes_old_SOURCES = cgstokes-old.cc
cgstokes_old_CPPFLAGS =				\
	$(AM_CPPFLAGS)				\
	$(SUPERLU_CPPFLAGS)			\
	$(UG_CPPFLAGS)				\
	$(ALUGRID_CPPFLAGS)
cgstokes_old_LDFLAGS =				\
	$(AM_LDFLAGS)				\
	$(SUPERLU_LDFLAGS)			\
	$(UG_LDFLAGS)				\
	$(ALUGRID_LDFLAGS)
cgstokes_old_LDADD =				\
	$(SUPERLU_LIBS)				\
	$(ALUGRID_LIBS)				\
	$(UG_LIBS)

endif

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
