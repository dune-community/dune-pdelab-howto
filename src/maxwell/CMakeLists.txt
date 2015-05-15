add_executable(modelproblem modelproblem.cc)
add_dune_alberta_flags(modelproblem)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vtk)
