# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
dist_noinst_DATA = poisson laplacedirichletccfv

noinst_PROGRAMS = poisson laplacedirichletccfv

poisson_SOURCES = poisson.cc
poisson_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
poisson_LDADD = $(MPI_LDFLAGS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
poisson_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
poisson_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

laplacedirichletccfv_SOURCES = laplacedirichletccfv.cc
laplacedirichletccfv_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
laplacedirichletccfv_LDADD = $(MPI_LDFLAGS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
laplacedirichletccfv_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
laplacedirichletccfv_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
