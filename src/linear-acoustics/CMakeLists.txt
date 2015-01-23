add_executable(heterogeneoussquare heterogeneoussquare.cc)
add_dune_alberta_flags(heterogeneoussquare)

set(CLEANFILES  dgfparser.log *.vtu *.png)
