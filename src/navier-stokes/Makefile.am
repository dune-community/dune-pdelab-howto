# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
dist_noinst_DATA =				\
	cgstokes.ini				\
	cgstokes_instat.ini			\
	grids/lshape.msh			\
	grids/pipe.msh				\
	grids/turbtube.msh			\
	grids/turbtube2d.msh

GRIDDIM=2

if SUPERLU
SUPERLUPROGRAMS = dgstokes \
	cgstokes \
	cgstokes_instat
endif

noinst_PROGRAMS = $(SUPERLUPROGRAMS)

noinst_HEADERS =				\
	cgstokes_initial.hh			\
	sproblemA.hh

#	novlp_test

#	dnaplfvamg

# novlp_test_SOURCES = novlp_test.cc
# novlp_test_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) \
# $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) \
# $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(PARMETIS_CPPFLAGS)
# novlp_test_LDADD = $(PARMETIS_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) \
# $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
# $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) \
# $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
# novlp_test_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS) $(PARMETIS_LDFLAGS)
# novlp_test_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

if SUPERLU
cgstokes_SOURCES = cgstokes.cc
cgstokes_CXXFLAGS = $(AM_CXXFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) -DDUNE_DEVEL_MODE $(ALUGRID_CPPFLAGS)
cgstokes_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(UG_LDFLAGS) $(UG_LIBS)
cgstokes_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
cgstokes_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
endif

if SUPERLU
cgstokes_instat_SOURCES = cgstokes_instat.cc
cgstokes_instat_CXXFLAGS = $(AM_CXXFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) -DDUNE_DEVEL_MODE $(ALUGRID_CPPFLAGS)
cgstokes_instat_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(UG_LDFLAGS) $(UG_LIBS)
cgstokes_instat_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
cgstokes_instat_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
endif

if SUPERLU
dgstokes_SOURCES = dgstokes.cc
dgstokes_CXXFLAGS = $(AM_CXXFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) -DDUNE_DEVEL_MODE
dgstokes_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(UG_LDFLAGS) $(UG_LIBS)
dgstokes_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
dgstokes_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
endif


CLEANFILES = dgfparser.log q1*.vtu *.vtu *.png

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
