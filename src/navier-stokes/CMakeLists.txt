# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA                               
        cgstokes.ini                            
        cgstokes_instat.ini                     
        grids/lshape.msh                        
        grids/pipe.msh                          
        grids/turbtube.msh                      
        grids/turbtube2d.msh)

set(GRIDDIM 2)

if(SUPERLU_FOUND)
set(SUPERLUPROGRAMS  dgstokes 
        cgstokes 
        cgstokes_instat)
endif(SUPERLU_FOUND )

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS})

set(noinst_HEADERS                                 
        cgstokes_initial.hh                     
        sproblemA.hh)

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

if(SUPERLU_FOUND)
add_executable(cgstokes cgstokes.cc)

add_dune_ug_flags(cgstokes )
add_dune_mpi_flags(cgstokes )
target_link_libraries(cgstokes ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes )
add_dune_alugrid_flags(cgstokes )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${UG_LDFLAGS} ${UG_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(cgstokes ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes )
endif(SUPERLU_FOUND )

if(SUPERLU_FOUND)
add_executable(cgstokes_instat cgstokes_instat.cc)

add_dune_ug_flags(cgstokes_instat )
add_dune_mpi_flags(cgstokes_instat )
target_link_libraries(cgstokes_instat ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes_instat )
add_dune_alugrid_flags(cgstokes_instat )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_instat_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${UG_LDFLAGS} ${UG_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_instat_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(cgstokes_instat ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes_instat )
endif(SUPERLU_FOUND )

if(SUPERLU_FOUND)
add_executable(dgstokes dgstokes.cc)

add_dune_ug_flags(dgstokes )
add_dune_mpi_flags(dgstokes )
target_link_libraries(dgstokes ${DUNE_LIBS})
add_dune_superlu_flags(dgstokes )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgstokes_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS} ${UG_LDFLAGS} ${UG_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgstokes_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(dgstokes ${DUNE_LIBS})
add_dune_superlu_flags(dgstokes )
endif(SUPERLU_FOUND )


set(CLEANFILES  dgfparser.log q1*.vtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
