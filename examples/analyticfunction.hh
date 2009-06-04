#include<math.h>
#include<dune/pdelab/common/function.hh>

template<typename GV, typename RF>
class U : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  U<GV,RF> > {
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,U<GV,RF> > B;

  U (const GV& gv) : B(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	typename Traits::DomainType center(0.0);
	center[0] = -80.0;
	center[0] = 20.0;
	center[0] = 0.0;
	center -= x;
	y = exp((-1.0/(25*25))*center.two_norm2());
  }
};
