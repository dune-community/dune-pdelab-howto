#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include"q1localfiniteelement.hh"

template<class D, class R>
class Q1LocalFiniteElementMap
  : public Dune::PDELab::
      SimpleLocalFiniteElementMap< Q1LocalFiniteElement<D,R> >
{};
