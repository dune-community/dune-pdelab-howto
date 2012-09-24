# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)
set(dist_noinst_DATA )

set(GRIDDIM 2)


set(noinst_PROGRAMS                                
        example01                               
        example02                               
        example03                               
        example04                               
        example05                               
        example06                               
        example07)

#       novlp_test

#       dnaplfvamg
set(noinst_HEADERS                                 
        example01_main.hh                       
        example01a_Q1.hh                        
        example01a_Q2.hh                        
        example01a_RT.hh                        
        example01a_operator.hh                  
        example01b_Q2.hh                        
        example01b_operator.hh                  
        example02_Q1.hh                         
        example02_bcextension.hh                
        example02_bctype.hh                     
        example02_operator.hh                   
        example03_Q2.hh                         
        example03_bcextension.hh                
        example03_bctype.hh                     
        example03_operator.hh                   
        example03_toperator.hh                  
        example04.hh                            
        example04_operator.hh                   
        example05_Q1Q1.hh                       
        example05_Q2Q2.hh                       
        example05_initial.hh                    
        example05_operator.hh                   
        example05_toperator.hh                  
        example06_Q1Q1.hh                       
        example06_bctype.hh                     
        example07_adaptivity.hh)

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


add_executable(example01 example01.cc)

add_dune_alberta_flags(example01 )
add_dune_amiramesh_flags(example01 )
add_dune_ug_flags(example01 )
add_dune_mpi_flags(example01 )
target_link_libraries(example01 ${DUNE_LIBS})
add_dune_superlu_flags(example01 )
add_dune_alugrid_flags(example01 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example01_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example01_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example01 ${DUNE_LIBS})
add_dune_superlu_flags(example01 )

add_executable(example02 example02.cc)

add_dune_alberta_flags(example02 )
add_dune_amiramesh_flags(example02 )
add_dune_ug_flags(example02 )
add_dune_mpi_flags(example02 )
target_link_libraries(example02 ${DUNE_LIBS})
add_dune_superlu_flags(example02 )
add_dune_alugrid_flags(example02 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example02_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example02_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example02 ${DUNE_LIBS})
add_dune_superlu_flags(example02 )

add_executable(example03 example03.cc)

add_dune_alberta_flags(example03 )
add_dune_amiramesh_flags(example03 )
add_dune_ug_flags(example03 )
add_dune_mpi_flags(example03 )
target_link_libraries(example03 ${DUNE_LIBS})
add_dune_superlu_flags(example03 )
add_dune_alugrid_flags(example03 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example03_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example03_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example03 ${DUNE_LIBS})
add_dune_superlu_flags(example03 )

add_executable(example04 example04.cc)

add_dune_alberta_flags(example04 )
add_dune_amiramesh_flags(example04 )
add_dune_ug_flags(example04 )
add_dune_mpi_flags(example04 )
target_link_libraries(example04 ${DUNE_LIBS})
add_dune_superlu_flags(example04 )
add_dune_alugrid_flags(example04 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example04_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example04_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example04 ${DUNE_LIBS})
add_dune_superlu_flags(example04 )

add_executable(example05 example05.cc)

add_dune_alberta_flags(example05 )
add_dune_amiramesh_flags(example05 )
add_dune_ug_flags(example05 )
add_dune_mpi_flags(example05 )
target_link_libraries(example05 ${DUNE_LIBS})
add_dune_superlu_flags(example05 )
add_dune_alugrid_flags(example05 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example05_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example05_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example05 ${DUNE_LIBS})
add_dune_superlu_flags(example05 )

add_executable(example06 example06.cc)

add_dune_alberta_flags(example06 )
add_dune_amiramesh_flags(example06 )
add_dune_ug_flags(example06 )
add_dune_mpi_flags(example06 )
target_link_libraries(example06 ${DUNE_LIBS})
add_dune_superlu_flags(example06 )
add_dune_alugrid_flags(example06 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example06_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example06_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(example06 ${DUNE_LIBS})
add_dune_superlu_flags(example06 )

add_executable(example07 example07.cc)
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example07_ LDADD                               
#        ${ALUGRID_LIBS}                         
#        ${UG_LIBS}                              
#        ${DUNEMPILIBS}                          
#        ${SUPERLU_LIBS}                         
#        ${LDADD})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(example07_ LDFLAGS ${AM_LDFLAGS}               
#        ${SUPERLU_LDFLAGS}                      
#        ${DUNEMPILDFLAGS}                       
#        ${UG_LDFLAGS}                           
#        ${ALUGRID_LDFLAGS})

add_dune_ug_flags(example07 )
add_dune_mpi_flags(example07 )
target_link_libraries(example07 ${DUNE_LIBS})
add_dune_superlu_flags(example07 )
add_dune_alugrid_flags(example07 )

set(CLEANFILES  dgfparser.log q1*.vtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
