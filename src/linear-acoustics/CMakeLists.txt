# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA                               
        splitdomain.geo                         
        splitdomain.msh)

set(GRIDDIM 2)

if(SUPERLU_FOUND)
set(SUPERLUPROGRAMS )
endif(SUPERLU_FOUND )

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS} 
        heterogeneoussquare)

if(SUPERLU_FOUND)
endif(SUPERLU_FOUND )

add_executable(heterogeneoussquare heterogeneoussquare.cc)

add_dune_alberta_flags(heterogeneoussquare )
add_dune_amiramesh_flags(heterogeneoussquare )
add_dune_ug_flags(heterogeneoussquare )
add_dune_mpi_flags(heterogeneoussquare )
target_link_libraries(heterogeneoussquare ${DUNE_LIBS})
add_dune_superlu_flags(heterogeneoussquare )
add_dune_alugrid_flags(heterogeneoussquare )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(heterogeneoussquare_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(heterogeneoussquare_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(heterogeneoussquare ${DUNE_LIBS})
add_dune_superlu_flags(heterogeneoussquare )

set(CLEANFILES  dgfparser.log *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
