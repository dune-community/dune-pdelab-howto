#include<dune/geometry/type.hh>
#include<dune/localfunctions/common/localfiniteelementtraits.hh>

#include"q1localbasis.hh"
#include"q1localcoefficients.hh"
#include"q1localinterpolation.hh"

template<class D, class R>
class Q1LocalFiniteElement
{
  Q1LocalBasis<D,R> basis;
  Q1LocalCoefficients coefficients;
  Q1LocalInterpolation<Q1LocalBasis<D,R> > interpolation;
  Dune::GeometryType gt;
public:
  typedef Dune::LocalFiniteElementTraits<Q1LocalBasis<D,R>,
				Q1LocalCoefficients,
	   			Q1LocalInterpolation<Q1LocalBasis<D,R> > > Traits;

  Q1LocalFiniteElement () { gt.makeQuadrilateral(); }

  const typename Traits::LocalBasisType& localBasis () const 
  {
	return basis;
  }
  
  const typename Traits::LocalCoefficientsType& localCoefficients () const 
  {
	return coefficients;
  }
  
  const typename Traits::LocalInterpolationType& localInterpolation () const 
  {
	return interpolation;
  }
	
  Dune::GeometryType type () const { return gt; }

  Q1LocalFiniteElement* clone () const {
    return new Q1LocalFiniteElement(*this);
  }
};
