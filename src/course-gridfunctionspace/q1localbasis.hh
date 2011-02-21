#include<dune/common/fvector.hh>
#include<dune/localfunctions/common/localbasis.hh>
template<class D, class R>
class Q1LocalBasis
{
public:
  typedef Dune::LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,
    Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,2>, 1> Traits;

  unsigned int size () const { return 4; }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (
				 const typename Traits::DomainType& in,
		   std::vector<typename Traits::RangeType>& out) const { 
	out.resize(4);
	out[0] = (1-in[0])*(1-in[1]); out[1] = (  in[0])*(1-in[1]);
	out[2] = (1-in[0])*(  in[1]); out[3] = (  in[0])*(  in[1]);
  }
  
  //! \brief Evaluate Jacobian of all shape functions
  inline void 
  evaluateJacobian (const typename Traits::DomainType& in,
		   std::vector<typename Traits::JacobianType>& out) const {  
	out.resize(4);
	out[0][0][0] = in[1]-1; out[0][0][1] = in[0]-1; /*@\label{q1b:grad0}@*/
	out[1][0][0] = 1-in[1]; out[1][0][1] = -in[0]; 
	out[2][0][0] =  -in[1]; out[2][0][1] = 1-in[0]; 
	out[3][0][0] =   in[1]; out[3][0][1] = in[0];   /*@\label{q1b:grad3}@*/
  }
  
  //! \brief Polynomial order of the shape functions
  unsigned int order () const  {
	return 1;
  }
};
