template<class GridView>
void do_something(const GridView &grid)
{
    // iterate over the leaf
    for (auto iterator = grid.template begin<0>(); iterator != grid.template end<0>(); ++iterator)
    {
      auto &entity = *iterator;
    }
}
