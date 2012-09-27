set(utilitydir  ${CMAKE_INSTALL_INCLUDEDIR}/src/utility)
set(utility_HEADERS  
        permeability_generator.hh 
        gridexamples.hh 
        basicunitcube.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${utility_HEADERS} DESTINATION ${utilitydir})