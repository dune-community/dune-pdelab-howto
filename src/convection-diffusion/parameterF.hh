#ifndef DUNE_PARAMETERF_HH
#define DUNE_PARAMETERF_HH

#include<math.h>

static char DurlofskyField[401] = 
"\
X..XX....XX....X.X..\
....X...XXX.........\
.....X..............\
.........X.....X...X\
......X....X......X.\
...X...X.....X......\
........X....XX..X..\
.................XX.\
X.X....XXX..X.......\
.........X.X.....X..\
X.........X....X....\
....X.XX.X.........X\
..X..XXX..X...X...X.\
X.X..X.XX...........\
.........X.....XX...\
.XX..XX....X......X.\
.XX....X......X...X.\
....X.............X.\
.........XX..X...X.X\
...........X.X.X....\
";

template<typename GV, typename RF>
class ParameterF
{
private:
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;


  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
	RF k;
	int X,Y,N;
	
	X = (int) (xglobal[0]*20); X = std::max(X,0); X = std::min(X,19);
	Y = (int) (xglobal[1]*20); Y = std::max(Y,0); Y = std::min(Y,19);
	N = (19-Y)*20+X;

	if (DurlofskyField[N]=='X')
	  k = 1E-6;
	else
	  k = 1.0;

    typename Traits::PermTensorType I;
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
	return 1.0; 
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1E-6 || xglobal[0]>1.0-1E-6)
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
	else
	  return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if (xglobal[0]<1E-6 )
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


#endif // DUNE_PARAMETERF_HH
