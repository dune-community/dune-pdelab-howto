auto it = gridview.template begin<0>();

// iterate over intersection of current entity
auto iitend = gridview.iend(*it);
for (auto iit = gridview.ibegin(*it); iit != iitend; ++iit)
{
  // neighbor intersection
  if (iit->neighbor())
  {
      // do something ...
  }
  // boundary intersection
  if (iit->boundary())
  {
      // do something else ...
  }
}
