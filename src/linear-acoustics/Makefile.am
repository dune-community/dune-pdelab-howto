# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
dist_noinst_DATA =				\
	splitdomain.geo				\
	splitdomain.msh

GRIDDIM=2

if SUPERLU
SUPERLUPROGRAMS =
endif

noinst_PROGRAMS = $(SUPERLUPROGRAMS) \
	heterogeneoussquare

if SUPERLU
endif

heterogeneoussquare_SOURCES = heterogeneoussquare.cc
heterogeneoussquare_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
heterogeneoussquare_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
heterogeneoussquare_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
heterogeneoussquare_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

CLEANFILES = dgfparser.log *.vtu *.png

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
