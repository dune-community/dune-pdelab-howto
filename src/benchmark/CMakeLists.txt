
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
add_executable(cgstokes_old EXCLUDE_FROM_ALL cgstokes-old.cc)

add_dune_ug_flags(cgstokes_old )
target_link_libraries(cgstokes_old ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes_old )
add_dune_alugrid_flags(cgstokes_old )

list(APPEND SUPERLUPROGRAMS cgstokes-benchmark)
add_executable(cgstokes-benchmark EXCLUDE_FROM_ALL cgstokes.cc)
list(APPEND SUPERLUPROGRAMS cgstokes-benchmark)
add_dune_ug_flags(cgstokes-benchmark)
target_link_libraries(cgstokes-benchmark ${DUNE_LIBS})
add_dune_superlu_flags(cgstokes-benchmark)
add_dune_alugrid_flags(cgstokes-benchmark)


# Generate list of macors for each executable
set(velocityblocked_macros CGSTOKES_BLOCK_LEVEL=1)
set(allblocked_macros CGSTOKES_BLOCK_LEVEL=2)
set(macroblocks_macros DCGSTOKES_MACROBLOCKS)
set(velocityblocked-macroblocks_macros
  CGSTOKES_BLOCK_LEVEL=1
  CGSTOKES_MACROBLOCKS)
set(suffixes velocityblocked allblocked macroblocks velocityblocked-macroblocks)

foreach(suffix ${suffixes})
  list(APPEND SUPERLUPROGRAMS cgstokes-${suffix})
  add_executable(cgstokes-${suffix} EXCLUDE_FROM_ALL cgstokes.cc)
  set_property(TARGET cgstokes-${suffix} APPEND
    PROPERTY COMPILE_DEFINITIONS ${${suffix}_macros})
  add_dune_ug_flags(cgstokes-${suffix})
  target_link_libraries(cgstokes-${suffix} ${DUNE_LIBS})
  add_dune_superlu_flags(cgstokes-${suffix})
  add_dune_alugrid_flags(cgstokes-${suffix})
endforeach(suffix ${suffixes})
endif(SUPERLU_FOUND )

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
