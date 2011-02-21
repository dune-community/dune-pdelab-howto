// interface: <dune/localfunctions/common/interpolation.hh>
template<class LB>
class Q1LocalInterpolation
{
public:

  //! \brief Local interpolation of a function
  template<typename F, typename C>
  void interpolate (const F& f, std::vector<C>& out) const {
	typename LB::Traits::DomainType x;
	typename LB::Traits::RangeType y;
	
	out.resize(4);
	x[0] = 0.0; x[1] = 0.0; f.evaluate(x,y); out[0] = y;
	x[0] = 1.0; x[1] = 0.0; f.evaluate(x,y); out[1] = y;
	x[0] = 0.0; x[1] = 1.0; f.evaluate(x,y); out[2] = y;
	x[0] = 1.0; x[1] = 1.0; f.evaluate(x,y); out[3] = y;
  }
};
