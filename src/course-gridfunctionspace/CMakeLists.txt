# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA )

set(GRIDDIM 2)

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS} 
        parallelgfs 
        q1gridfunctionspacemain 
        q1interpolatemain 
        q1interpolationerrormain 
        q2interpolationerrormain 
        q1constrainedinterpolatemain 
        thinterpolatemain)

set(noinst_HEADERS                                 
        analyticfunction.hh                     
        l2interpolationerror.hh                 
        q1constrainedinterpolate.hh             
        q1constraints.hh                        
        q1constraintsparameters.hh              
        q1gridfunctionspace.hh                  
        q1interpolate.hh                        
        q1interpolationerror.hh                 
        q1localbasis.hh                         
        q1localcoefficients.hh                  
        q1localfiniteelement.hh                 
        q1localfiniteelementmap.hh              
        q1localinterpolation.hh                 
        q2interpolationerror.hh                 
        thinterpolate.hh                        
        thvelocity.hh)

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

add_executable(parallelgfs parallelgfs.cc)

add_dune_alberta_flags(parallelgfs )
add_dune_amiramesh_flags(parallelgfs )
add_dune_ug_flags(parallelgfs )
add_dune_mpi_flags(parallelgfs )
target_link_libraries(parallelgfs ${DUNE_LIBS})
add_dune_superlu_flags(parallelgfs )
add_dune_alugrid_flags(parallelgfs )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(parallelgfs_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(parallelgfs_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(parallelgfs ${DUNE_LIBS})
add_dune_superlu_flags(parallelgfs )

add_executable(q1gridfunctionspacemain q1gridfunctionspacemain.cc)

add_dune_alberta_flags(q1gridfunctionspacemain )
add_dune_amiramesh_flags(q1gridfunctionspacemain )
add_dune_ug_flags(q1gridfunctionspacemain )
add_dune_mpi_flags(q1gridfunctionspacemain )
add_dune_alugrid_flags(q1gridfunctionspacemain )
target_link_libraries(q1gridfunctionspacemain ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1gridfunctionspacemain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})

add_executable(q1interpolatemain q1interpolatemain.cc)

add_dune_alberta_flags(q1interpolatemain )
add_dune_amiramesh_flags(q1interpolatemain )
add_dune_ug_flags(q1interpolatemain )
add_dune_mpi_flags(q1interpolatemain )
target_link_libraries(q1interpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(q1interpolatemain )
add_dune_alugrid_flags(q1interpolatemain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1interpolatemain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1interpolatemain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(q1interpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(q1interpolatemain )

add_executable(q1interpolationerrormain q1interpolationerrormain.cc)

add_dune_alberta_flags(q1interpolationerrormain )
add_dune_amiramesh_flags(q1interpolationerrormain )
add_dune_ug_flags(q1interpolationerrormain )
add_dune_mpi_flags(q1interpolationerrormain )
target_link_libraries(q1interpolationerrormain ${DUNE_LIBS})
add_dune_superlu_flags(q1interpolationerrormain )
add_dune_alugrid_flags(q1interpolationerrormain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1interpolationerrormain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1interpolationerrormain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(q1interpolationerrormain ${DUNE_LIBS})
add_dune_superlu_flags(q1interpolationerrormain )

add_executable(q2interpolationerrormain q2interpolationerrormain.cc)

add_dune_alberta_flags(q2interpolationerrormain )
add_dune_amiramesh_flags(q2interpolationerrormain )
add_dune_ug_flags(q2interpolationerrormain )
add_dune_mpi_flags(q2interpolationerrormain )
target_link_libraries(q2interpolationerrormain ${DUNE_LIBS})
add_dune_superlu_flags(q2interpolationerrormain )
add_dune_alugrid_flags(q2interpolationerrormain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q2interpolationerrormain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q2interpolationerrormain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(q2interpolationerrormain ${DUNE_LIBS})
add_dune_superlu_flags(q2interpolationerrormain )

add_executable(q1constrainedinterpolatemain q1constrainedinterpolatemain.cc)

add_dune_alberta_flags(q1constrainedinterpolatemain )
add_dune_amiramesh_flags(q1constrainedinterpolatemain )
add_dune_ug_flags(q1constrainedinterpolatemain )
add_dune_mpi_flags(q1constrainedinterpolatemain )
target_link_libraries(q1constrainedinterpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(q1constrainedinterpolatemain )
add_dune_alugrid_flags(q1constrainedinterpolatemain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1constrainedinterpolatemain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(q1constrainedinterpolatemain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(q1constrainedinterpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(q1constrainedinterpolatemain )

add_executable(thinterpolatemain thinterpolatemain.cc)

add_dune_alberta_flags(thinterpolatemain )
add_dune_amiramesh_flags(thinterpolatemain )
add_dune_ug_flags(thinterpolatemain )
add_dune_mpi_flags(thinterpolatemain )
target_link_libraries(thinterpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(thinterpolatemain )
add_dune_alugrid_flags(thinterpolatemain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(thinterpolatemain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(thinterpolatemain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(thinterpolatemain ${DUNE_LIBS})
add_dune_superlu_flags(thinterpolatemain )

set(CLEANFILES  dgfparser.log q1*.vtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
