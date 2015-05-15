add_executable(heterogeneoussquare heterogeneoussquare.cc)
add_dune_alberta_flags(heterogeneoussquare)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vtk)

set(CLEANFILES  dgfparser.log *.vtu *.png)
