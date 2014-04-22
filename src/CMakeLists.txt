set(SUBDIRS
	boilerplatetutorial
        convection-diffusion
        course-examples
        course-gridfunctionspace
        linear-acoustics
        maxwell
        navier-stokes
        two-phase
        utility)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
