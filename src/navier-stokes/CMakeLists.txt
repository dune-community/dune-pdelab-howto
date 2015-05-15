if(SUPERLU_FOUND)
add_executable(cgstokes cgstokes.cc)
add_executable(cgstokes_instat cgstokes_instat.cc)
add_executable(dgstokes dgstokes.cc)
endif()

dune_symlink_to_source_files(FILES cgstokes.ini cgstokes_instat.ini grids)
