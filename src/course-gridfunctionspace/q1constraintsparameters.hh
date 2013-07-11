#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

class BParam
  : public Dune::PDELab::DirichletConstraintsParameters                             /*@\label{bcp:base}@*/
{
public:
  template<typename I>
  bool isDirichlet(const I & intersection,                                          /*@\label{bcp:name}@*/
    const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global(coord);

    if (xg[0]<1E-6 && xg[1]>0.25 && xg[1]<0.75 )
      return false; // no Dirichlet
    return true;  // Dirichlet
  }
};
