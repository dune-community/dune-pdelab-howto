#include <dune/grid/somegrid.hh>

void iterate_the_grid()
{
    // we have a grid
    Dune::SomeGrid grid(parameters);

    // iterate over level 2
    for (auto it = grid.levelGridView(2).template begin<0>(); it != grid.levelGridView(2).template end<0>(); ++it) {
        ...
    }
    // iterate over the leaf
    for (auto it = grid.leafGridView().template begin<0>(); it != grid.leafGridView().template end(); ++it) {
        ...
    }
}
