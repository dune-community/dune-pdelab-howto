# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA )

set(GRIDDIM 2)

if(SUPERLU_FOUND)
set(SUPERLUPROGRAMS  dnaplmfd)
endif(SUPERLU_FOUND )

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS}
        dnaplfv)

set(noinst_HEADERS  twophasemfd.hh)

#       novlp_test

#       dnaplfvamg

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

add_executable(dnaplfv dnaplfv.cc)

add_dune_alberta_flags(dnaplfv )
add_dune_amiramesh_flags(dnaplfv )
add_dune_parmetis_flags(dnaplfv )
add_dune_ug_flags(dnaplfv )
add_dune_mpi_flags(dnaplfv )
target_link_libraries(dnaplfv ${DUNE_LIBS})
add_dune_superlu_flags(dnaplfv )
add_dune_alugrid_flags(dnaplfv )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dnaplfv_ LDADD ${PARMETIS_LIBS} ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dnaplfv_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS} ${PARMETIS_LDFLAGS})

target_link_libraries(dnaplfv ${DUNE_LIBS})
add_dune_superlu_flags(dnaplfv )

# dnaplmfd is broken (there is no IntersectionIndexSet for the new backend structure)
# if(SUPERLU_FOUND)
# add_executable(dnaplmfd dnaplmfd.cc)
#
# add_dune_alberta_flags(dnaplmfd )
# add_dune_amiramesh_flags(dnaplmfd )
# add_dune_ug_flags(dnaplmfd )
# add_dune_mpi_flags(dnaplmfd )
# target_link_libraries(dnaplmfd ${DUNE_LIBS})
# add_dune_superlu_flags(dnaplmfd )
# add_dune_alugrid_flags(dnaplmfd )
# # The following lines are comment out as the translation is not clear. TODO review
# #
# # set(dnaplmfd_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# # The following lines are comment out as the translation is not clear. TODO review
# #
# # set(dnaplmfd_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})
#
# target_link_libraries(dnaplmfd ${DUNE_LIBS})
# add_dune_superlu_flags(dnaplmfd )
# endif(SUPERLU_FOUND )


set(CLEANFILES  dgfparser.log q1*.vtu *.pvtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
