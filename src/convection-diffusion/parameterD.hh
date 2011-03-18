#ifndef DUNE_PARAMETERD_HH
#define DUNE_PARAMETERD_HH

#include<math.h>
#include"../utility/permeability_generator.hh"

template<typename GV, typename RF>
class ParameterD
{

private:

  const GV& gv;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;

  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;


  ParameterD(
			 const GV& gv_,
			 Dune::FieldVector<double,GV::dimension> correlation_length,
			 double variance = 1.0, 
			 double mean = 0.0, 
			 long modes = 1000, 
			 long seed = -1083
			 )
  :
	gv(gv_), 
	is(gv.indexSet()), 
	perm(is.size(0))
  {
	typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	typedef typename Traits::DomainFieldType DF;
	const int dim = GV::dimension;
	double mink=1E100;
	double maxk=-1E100;
	
	
	EberhardPermeabilityGenerator<GV::dimension> field(correlation_length,variance,mean,modes,seed);
	
	for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
	  {
		int id = is.index(*it);
        Dune::GeometryType gt = it->geometry().type();
        Dune::FieldVector<DF,dim> localcenter =
          Dune::GenericReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
		perm[id]=field.eval(globalcenter);
		mink = std::min(mink,log10(perm[id]));
		maxk = std::max(maxk,log10(perm[id]));
	  }
	std::cout << "log10(mink)=" << mink << " log10(maxk)=" << maxk << std::endl;
  }


  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
	RF k = 0;
	k = perm[is.index(e)];

    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          I[i][i] = k;
        else
          I[i][j] = 0.0;
	
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType 
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
	return 0.0; 
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
	if( xglobal[0] < 1E-6 || xglobal[0] > 1.0-1E-6 )
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
	else
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
	if( xglobal[0] < 1E-6 )
	  return 1.0;
	else
	  return 0.0;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType 
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};


#endif // DUNE_PARAMETERD_HH
