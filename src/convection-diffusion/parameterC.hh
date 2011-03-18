#ifndef DUNE_PARAMETERC_HH
#define DUNE_PARAMETERC_HH

template<typename GV, typename RF>
class ParameterC
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


  ParameterC()
  {
	K000=20.0;
	K001=0.002;
	K010=0.2;
	K011=2000.0;
	K100=1000.0;
	K101=0.001;
	K110=0.1;
	K111=10.0;
	width=1.0/8.0;
  }


  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    typename Traits::DomainType xglobal = e.geometry().global(x);

    int ix,iy,iz;

    ix=((int)floor(xglobal[0]/width))%2;
    iy=((int)floor(xglobal[1]/width))%2;
	if (GV::dimension>2)
	  iz=((int)floor(xglobal[2]/width))%2;
	else
	  iz=0;

	RF k = 0;
    if ( iz==0 && iy==0 && ix==0 ) k=K000;
    if ( iz==0 && iy==0 && ix==1 ) k=K001;
    if ( iz==0 && iy==1 && ix==0 ) k=K010;
    if ( iz==0 && iy==1 && ix==1 ) k=K011;
    if ( iz==1 && iy==0 && ix==0 ) k=K100;
    if ( iz==1 && iy==0 && ix==1 ) k=K101;
    if ( iz==1 && iy==1 && ix==0 ) k=K110;
    if ( iz==1 && iy==1 && ix==1 ) k=K111;

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


#endif // DUNE_PARAMETERC_HH
