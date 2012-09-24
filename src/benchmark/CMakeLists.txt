
set(MOSTLYCLEANFILES )
set(SUPERLUPROGRAMS )

set(noinst_PROGRAMS  ${SUPERLUPROGRAMS})

set(AM_CPPFLAGS                                    
        ${AM_CPPFLAGS}                           
        ${DUNE_CPPFLAGS}                        
        ${DUNEMPIFLAGS})

set(AM_LDFLAGS                                     
        ${AM_LDFLAGS}                            
        ${DUNE_LDFLAGS}                         
        ${DUNEMPILDFLAGS})

set(LDADD                                          
        ${DUNE_LIBS}                            
        ${DUNEMPILIBS}                          
        ${LDADD})

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# TRUNK TARGETS FOLLOW

if(SUPERLU_FOUND)

list(APPEND MOSTLYCLEANFILES cgstokes_*.vtu cgstokes_*_timings.txt)

list(APPEND SUPERLUPROGRAMS cgstokes-old)
add_executable(cgstokes_old cgstokes-old.cc)

add_dune_ug_flags(cgstokes_old )
target_link_libraries(cgstokes_old ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes_old )
add_dune_alugrid_flags(cgstokes_old )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_old_ LDFLAGS                          
#        ${AM_LDFLAGS}                           
#        ${SUPERLU_LDFLAGS}                      
#        ${UG_LDFLAGS}                           
#        ${ALUGRID_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(cgstokes_old_ LDADD                            
#        ${SUPERLU_LIBS}                         
#        ${ALUGRID_LIBS}                         
#        ${UG_LIBS})

endif(SUPERLU_FOUND )

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
