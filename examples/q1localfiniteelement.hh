#include<dune/common/geometrytype.hh>
#include<dune/finiteelements/common/localfiniteelement.hh>

template<class D, class R>
class Q1LocalFiniteElement : public Dune::LocalFiniteElementInterface< 
  Dune::LocalFiniteElementTraits<Q1LocalBasis<D,R>,Q1LocalCoefficients,
								 Q1LocalInterpolation<Q1LocalBasis<D,R> > >, 
  Q1LocalFiniteElement<D,R> >
{
public:
  typedef Dune::LocalFiniteElementTraits<Q1LocalBasis<D,R>,
										 Q1LocalCoefficients,
	   				 Q1LocalInterpolation<Q1LocalBasis<D,R> > > Traits;

  Q1LocalFiniteElement () {
	gt.makeQuadrilateral();
  }

  const typename Traits::LocalBasisType& localBasis () const {
	return basis;
  }
  
  const typename Traits::LocalCoefficientsType& localCoefficients () const {
	return coefficients;
  }
  
  const typename Traits::LocalInterpolationType& localInterpolation () const {
	return interpolation;
  }
	
  Dune::GeometryType type () const {
	return gt;
  }

private:
  Q1LocalBasis<D,R> basis;
  Q1LocalCoefficients coefficients;
  Q1LocalInterpolation<Q1LocalBasis<D,R> > interpolation;
  Dune::GeometryType gt;
};
