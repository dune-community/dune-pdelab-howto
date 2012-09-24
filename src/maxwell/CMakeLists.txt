# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA                               
        mysphere.geo                            
        mysphere1.msh                           
        mysphere2.msh)

set(GRIDDIM 2)

if(SUPERLU_FOUND)
set(SUPERLUPROGRAMS )
endif(SUPERLU_FOUND )

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS} 
        modelproblem)

if(SUPERLU_FOUND)
endif(SUPERLU_FOUND )

add_executable(modelproblem modelproblem.cc)

add_dune_alberta_flags(modelproblem )
add_dune_amiramesh_flags(modelproblem )
add_dune_ug_flags(modelproblem )
add_dune_mpi_flags(modelproblem )
target_link_libraries(modelproblem ${DUNE_LIBS})
add_dune_superlu_flags(modelproblem )
add_dune_alugrid_flags(modelproblem )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(modelproblem_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(modelproblem_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(modelproblem ${DUNE_LIBS})
add_dune_superlu_flags(modelproblem )

set(CLEANFILES  dgfparser.log *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
