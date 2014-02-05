# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

SUBDIRS =

GRIDDIM=2

if EIGEN
EIGENPROGRAMS =	\
    laplacedirichletccfv-eigen
endif

if SUPERLU
SUPERLUPROGRAMS =				\
    dgdiffusion					\
	dgdiffusion-hp				\
	dgoverlapping				\
	rt0main
endif

noinst_PROGRAMS =				\
    $(EIGENPROGRAMS)				\
    $(SUPERLUPROGRAMS)				\
    $(ALUGRIDPROGRAMS)                          \
	overlappingsinglephaseflowccfv		\
	overlappingsinglephaseflow		\
	nonoverlappingsinglephaseflow		\
	nonlineardiffusion			\
	transporttest				\
	instationarytest			\
	poisson		\
	poisson_HN		\
	laplacedirichletccfv			\
	reentrantcorner				\
	diffusion				\
	scalabilitytest				\
	dgldomain				\
	ldomain

# mimetic is broken (probably irreparably...)


noinst_HEADERS =				\
	l2interpolationerror.hh			\
	parameterA.hh				\
	parameterB.hh				\
	parameterC.hh				\
	parameterD.hh				\
	parameterE.hh				\
	parameterF.hh				\
	problemA.hh				\
	problemB.hh				\
	problemC.hh				\
	problemD.hh				\
	problemE.hh				\
	problemF.hh				\
	reentrantcornerproblem.hh

dist_noinst_DATA =				\
	grids/2dgrid.al				\
	grids/2dreentrantcorner.dgf		\
	grids/2dsimplex.alu			\
	grids/3dgrid.al				\
	grids/3drefinedcube.alu			\
	grids/CuKav.geo				\
	grids/CuKav.msh				\
	grids/circle.geo			\
	grids/circle1storder.msh		\
	grids/circle2ndorder.msh		\
	grids/cube.geo				\
	grids/cube.hexa				\
	grids/cube.tetra			\
	grids/cube1045.msh			\
	grids/cube3205.msh			\
	grids/cube99.msh			\
	grids/curved2d.geo			\
	grids/curved2d.msh			\
	grids/ldomain.al			\
	grids/ldomain.geo			\
	grids/ldomain.msh			\
	grids/ldomain2.msh			\
	grids/pyramid.geo			\
	grids/pyramid.msh			\
	grids/pyramid1storder.msh		\
	grids/pyramid2ndorder.msh		\
	grids/pyramid4.msh			\
	grids/telescope.geo			\
	grids/telescope.msh			\
	grids/telescope1storder.msh		\
	grids/telescope2ndorder.msh

if SUPERLU
dgdiffusion_SOURCES = dgdiffusion.cc
dgdiffusion_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
dgdiffusion_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
dgdiffusion_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
dgdiffusion_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

dgdiffusion_hp_SOURCES = dgdiffusion-hp.cc
dgdiffusion_hp_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
dgdiffusion_hp_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
dgdiffusion_hp_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
dgdiffusion_hp_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)
endif

dgoverlapping_SOURCES = dgoverlapping.cc
dgoverlapping_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
dgoverlapping_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
dgoverlapping_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
dgoverlapping_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

overlappingsinglephaseflowccfv_SOURCES = overlappingsinglephaseflowccfv.cc
overlappingsinglephaseflowccfv_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
overlappingsinglephaseflowccfv_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
overlappingsinglephaseflowccfv_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
overlappingsinglephaseflowccfv_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

overlappingsinglephaseflow_SOURCES = overlappingsinglephaseflow.cc
overlappingsinglephaseflow_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
overlappingsinglephaseflow_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
overlappingsinglephaseflow_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
overlappingsinglephaseflow_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

nonoverlappingsinglephaseflow_SOURCES = nonoverlappingsinglephaseflow.cc
nonoverlappingsinglephaseflow_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
nonoverlappingsinglephaseflow_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
nonoverlappingsinglephaseflow_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
nonoverlappingsinglephaseflow_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

poisson_SOURCES = poisson.cc
poisson_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
poisson_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
poisson_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
poisson_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

poisson_HN_SOURCES = poisson.cc
poisson_HN_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS) -DHANGING_NODES_REFINEMENT
poisson_HN_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
poisson_HN_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
poisson_HN_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

nonlineardiffusion_SOURCES = nonlineardiffusion.cc
nonlineardiffusion_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
nonlineardiffusion_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
nonlineardiffusion_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
nonlineardiffusion_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

transporttest_SOURCES = transporttest.cc
transporttest_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
transporttest_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
transporttest_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
transporttest_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

instationarytest_SOURCES = instationarytest.cc
instationarytest_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
instationarytest_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
instationarytest_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
instationarytest_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

# mimetic_SOURCES = mimetic.cc
# mimetic_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
# mimetic_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
# mimetic_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
# mimetic_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

laplacedirichletccfv_SOURCES = laplacedirichletccfv.cc
laplacedirichletccfv_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS)
laplacedirichletccfv_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(UG_LDFLAGS) $(UG_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

if EIGEN
laplacedirichletccfv_eigen_SOURCES = laplacedirichletccfv.cc
laplacedirichletccfv_eigen_CPPFLAGS = $(laplacedirichletccfv_CPPFLAGS) -DUSE_EIGEN
laplacedirichletccfv_eigen_CXXFLAGS = $(EIGEN_CXXFLAGS)
laplacedirichletccfv_eigen_LDADD = $(laplacedirichletccfv_LDADD)
endif

reentrantcorner_SOURCES = reentrantcorner.cc
reentrantcorner_CXXFLAGS = $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
reentrantcorner_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
reentrantcorner_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
reentrantcorner_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

if SUPERLU
rt0main_SOURCES = rt0main.cc
rt0main_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
rt0main_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
rt0main_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
endif

diffusion_SOURCES = diffusion.cc
diffusion_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
diffusion_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
diffusion_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
diffusion_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

scalabilitytest_SOURCES = scalabilitytest.cc
scalabilitytest_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
scalabilitytest_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
scalabilitytest_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
scalabilitytest_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

ldomain_SOURCES = ldomain.cc
ldomain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
ldomain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
ldomain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
ldomain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

dgldomain_SOURCES = dgldomain.cc
dgldomain_CXXFLAGS = $(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(ALUGRID_CPPFLAGS)
dgldomain_LDADD = $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
dgldomain_LDFLAGS= $(AM_LDFLAGS) $(SUPERLU_LDFLAGS)
dgldomain_CPPFLAGS=$(AM_CPPFLAGS) $(SUPERLU_CPPFLAGS)

CLEANFILES = poisson_*.vtu reentrant_*.vtu dgfparser.log q1*.vtu *.vtu *.png

#============================================================
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
