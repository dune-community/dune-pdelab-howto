#include<dune/pdelab/common/function.hh>
template<typename GV, typename R>
class V : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,R,GV::dimension>,
  V<GV,R> > {
public:
  typedef Dune::PDELab::
    AnalyticGridFunctionTraits<GV,R,GV::dimension> Traits;
  typedef Dune::PDELab::
    AnalyticGridFunctionBase<Traits,V<GV,R> > BaseT;

  V (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const {  
	y[0]=x[0]; for (int i=1; i<GV::dimension; i++) y[i]=x[i]*y[i-1];
  }
};
