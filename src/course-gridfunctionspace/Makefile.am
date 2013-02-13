# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
dist_noinst_DATA =

GRIDDIM=2

noinst_PROGRAMS = $(SUPERLUPROGRAMS) \
	parallelgfs \
	q1gridfunctionspacemain \
	q1interpolatemain \
	q1interpolationerrormain \
	q2interpolationerrormain \
	q1constrainedinterpolatemain \
	thinterpolatemain

noinst_HEADERS =				\
	analyticfunction.hh			\
	l2interpolationerror.hh			\
	q1constrainedinterpolate.hh		\
	q1constraints.hh			\
	q1constraintsparameters.hh		\
	q1gridfunctionspace.hh			\
	q1interpolate.hh			\
	q1interpolationerror.hh			\
	q1localbasis.hh				\
	q1localcoefficients.hh			\
	q1localfiniteelement.hh			\
	q1localfiniteelementmap.hh		\
	q1localinterpolation.hh			\
	q2interpolationerror.hh			\
	thinterpolate.hh			\
	thvelocity.hh

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

parallelgfs_SOURCES = parallelgfs.cc
parallelgfs_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
parallelgfs_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
parallelgfs_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
parallelgfs_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

q1gridfunctionspacemain_SOURCES = q1gridfunctionspacemain.cc
q1gridfunctionspacemain_CXXFLAGS = $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
q1gridfunctionspacemain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

q1interpolatemain_SOURCES = q1interpolatemain.cc
q1interpolatemain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
q1interpolatemain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
q1interpolatemain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
q1interpolatemain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

q1interpolationerrormain_SOURCES = q1interpolationerrormain.cc
q1interpolationerrormain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
q1interpolationerrormain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
q1interpolationerrormain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
q1interpolationerrormain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

q2interpolationerrormain_SOURCES = q2interpolationerrormain.cc
q2interpolationerrormain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
q2interpolationerrormain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
q2interpolationerrormain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
q2interpolationerrormain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

q1constrainedinterpolatemain_SOURCES = q1constrainedinterpolatemain.cc
q1constrainedinterpolatemain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
q1constrainedinterpolatemain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
q1constrainedinterpolatemain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
q1constrainedinterpolatemain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

thinterpolatemain_SOURCES = thinterpolatemain.cc
thinterpolatemain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
thinterpolatemain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
thinterpolatemain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
thinterpolatemain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

CLEANFILES = dgfparser.log q1*.vtu *.vtu *.png

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
