#ifndef DUNE_PARAMETERB_HH
#define DUNE_PARAMETERB_HH

template<typename GV, typename RF>
class ParameterB
{
  RF K000;
  RF K001;
  RF K010;
  RF K011;
  RF K100;
  RF K101;
  RF K110;
  RF K111;
  RF width;

  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;


  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          I[i][i] = 1.0;
        else
          I[i][j] = 0.0;

    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    typename Traits::DomainType xglobal = e.geometry().global(x);
	
    if( xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6)
      {
        v[0] = 0; 
		v[1] = 0;
      }
	
    if (xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6)
      {
        v[0] = -5.0; 
		v[1] = 0.0;
      }

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
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if (xglobal[0]>0.25 && xglobal[0]<0.375 && xglobal[1]>0.25 && xglobal[1]<0.375)
      return 50.0;
    else
      return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);

    if (xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6)
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;

    if (xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6)
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;

	return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= xglobal;
	return exp(-center.two_norm2());
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if( xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6 )
      {
        return 0;
      }
    if( xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6 )
      {
        return -5.0;
      }
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType 
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};


#endif // DUNE_PARAMETERB_HH
