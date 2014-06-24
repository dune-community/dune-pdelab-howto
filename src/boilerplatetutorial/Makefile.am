# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $
SUBDIRS =

noinst_PROGRAMS =				\
poisson_uniform					\
	poisson_adaptive			\
	advection_stationary			\
	advection_instationary			\
	heat_instationary


noinst_HEADERS =

dist_noinst_DATA =

poisson_uniform_SOURCES = poisson_uniform.cc
poisson_uniform_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
poisson_uniform_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
poisson_uniform_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
poisson_uniform_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

poisson_adaptive_SOURCES = poisson_adaptive.cc
poisson_adaptive_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
poisson_adaptive_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
poisson_adaptive_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
poisson_adaptive_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

advection_stationary_SOURCES = advection_stationary.cc
advection_stationary_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
advection_stationary_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
advection_stationary_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
advection_stationary_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

advection_instationary_SOURCES = advection_instationary.cc
advection_instationary_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
advection_instationary_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
advection_instationary_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
advection_instationary_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

heat_instationary_SOURCES = heat_instationary.cc
heat_instationary_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
heat_instationary_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
heat_instationary_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
heat_instationary_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

CLEANFILES = poisson_*.vtu reentrant_*.vtu dgfparser.log q1*.vtu *.vtu *.png

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules
