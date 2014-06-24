set(BOILERPLATEPROGRAMS
  poisson_uniform
  poisson_adaptive
  advection_stationary
  advection_instationary
  heat_instationary)

foreach(prog ${BOILERPLATEPROGRAMS})
  add_executable(${prog} ${prog}.cc)
  add_dune_superlu_flags(${prog})
  add_dune_alugrid_flags(${prog})
  add_dune_alberta_flags(${prog})
  add_dune_ug_flags(${prog})
  target_link_libraries(${prog} ${DUNE_LIBS})
endforeach(prog ${BOILERPLATEPROGRAMS})
