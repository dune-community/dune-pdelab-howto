#include<dune/common/fvector.hh>
#include<dune/finiteelements/common/localbasis.hh>

template<class D, class R>
class Q1LocalBasis : 
  public Dune::C1LocalBasisInterface<
  Dune::C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,
						   R,1,Dune::FieldVector<R,1>,
						   Dune::FieldVector<Dune::FieldVector<R,2>,1> >,
  Q1LocalBasis<D,R> >
{
public:
  typedef Dune::C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,
								   R,1,Dune::FieldVector<R,1>,
				  Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

  //! \brief number of shape functions
  unsigned int size () const
  {
	return 4;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& in,
						  std::vector<typename Traits::RangeType>& out) const
  { 
	out.resize(4);
	out[0] = (1-in[0])*(1-in[1]);
	out[1] = (  in[0])*(1-in[1]);
	out[2] = (1-in[0])*(  in[1]);
	out[3] = (  in[0])*(  in[1]);
  }
  
  //! \brief Evaluate Jacobian of all shape functions
  inline void 
  evaluateJacobian (const typename Traits::DomainType& in,
		   std::vector<typename Traits::JacobianType>& out) const
  {  
	out.resize(4);
	out[0][0][0] = in[1]-1; out[0][0][1] = in[0]-1; 
	out[1][0][0] = 1-in[1]; out[1][0][1] = -in[0]; 
	out[2][0][0] =  -in[1]; out[2][0][1] = 1-in[0]; 
	out[3][0][0] =   in[1]; out[3][0][1] = in[0]; 
  }
  
  //! \brief Polynomial order of the shape functions
  unsigned int order () const
  {
	return 1;
  }
};
