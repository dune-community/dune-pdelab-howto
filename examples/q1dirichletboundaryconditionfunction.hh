#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>

template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
  BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,B<GV> >
{

public:
  typedef Dune::PDELab::
  BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const I& ig, const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const {  
    y = 1; // Dirichlet
  }

  inline const GV& getGridView () {
    return gv;
  }
private:
  const GV& gv;
};
