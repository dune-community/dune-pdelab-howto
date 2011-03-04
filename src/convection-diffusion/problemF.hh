#ifndef DUNE_PARSOLVE_PROBLEMF_HH
#define DUNE_PARSOLVE_PROBLEMF_HH

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

// function for defining the diffusion tensor
template<typename GV, typename RF>
class k_F
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> >,
      k_F<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,k_F<GV,RF> > BaseT;

  k_F (const GV& gv_) : gv(gv_)
  { }

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  { 
	Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      xg = e.geometry().global(x);

	int X,Y,N;

	X = (int) (xg[0]*20); X = std::max(X,0); X = std::min(X,19);
	Y = (int) (xg[1]*20); Y = std::max(Y,0); Y = std::min(Y,19);
	N = (19-Y)*20+X;

	if (DurlofskyField[N]=='X')
	  y = 1E-6;
	else
	  y = 1.0;
  }
  
  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
  
private:
  const GV& gv;
};

// function for defining the diffusion tensor
template<typename GV, typename RF>
class K_F
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_F<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_F<GV,RF> > BaseT;

  K_F (const GV& gv_) : gv(gv_)
  {}
    
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  { 
	Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      xg = e.geometry().global(x);

	int X,Y,N;

	X = (int) (xg[0]*20); X = std::max(X,0); X = std::min(X,19);
	Y = (int) (xg[1]*20); Y = std::max(Y,0); Y = std::min(Y,19);
	N = (19-Y)*20+X;

	RF k;
	if (DurlofskyField[N]=='X')
	  k = 1E-6;
	else
	  k = 1.0;

    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
		  y[i][i] = k;
        else
          y[i][j] = 0.0;
  }
  
  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

private:
  const GV& gv;
};

// function for defining the source term
template<typename GV, typename RF>
class A0_F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  A0_F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,A0_F<GV,RF> > BaseT;

  A0_F (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class F_F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F_F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_F<GV,RF> > BaseT;

  F_F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0; 
  }
};


// constraints parameter class for selecting boundary condition type 
class BCTypeParam_F
  : public Dune::PDELab::FluxConstraintsParameters,
	public Dune::PDELab::DirichletConstraintsParameters
	/*@\label{bcp:base}@*/
{
public:

  template<typename I>
  bool isNeumann(
				   const I & intersection,   /*@\label{bcp:name}@*/
				   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
				   ) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
	  xg = intersection.geometry().global( coord );
	
    if (xg[0]<1E-6 || xg[0]>1.0-1E-6)
	  return false; // Dirichlet
	else
	  return true;
  }

  template<typename I>
  bool isDirichlet(
				   const I & intersection,   /*@\label{bcp:name}@*/
				   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
				   ) const
  {
	return !isNeumann( intersection, coord );
  }
};




// boundary grid function selecting boundary conditions 
template<typename GV>
class B_F
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<
                                                    GV,Dune::PDELab::DiffusionBoundaryCondition::Type,1,
                                                    Dune::FieldVector<
                                                      Dune::PDELab::DiffusionBoundaryCondition::Type,1> >,
                                                  B_F<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::DiffusionBoundaryCondition BC;
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,
    Dune::PDELab::DiffusionBoundaryCondition::Type,1,
    Dune::FieldVector<Dune::PDELab::DiffusionBoundaryCondition::Type,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_F<GV> > BaseT;

  B_F (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
	Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      xg = ig.geometry().global(x);

    if (xg[0]<1E-6 || xg[0]>1.0-1E-6)
      y = BC::Dirichlet;
	else
	  y = BC::Neumann;
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G_F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G_F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G_F<GV,RF> > BaseT;

  G_F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0;
    if (x[0]<1E-6)
	  y = 1;
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J_F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J_F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J_F<GV,RF> > BaseT;

  J_F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0;
	return;
  }
};

// flux as velocity field for the mixed method
template<typename GV, typename RF>
class V_F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
													  V_F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_F<GV,RF> > BaseT;

  V_F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
    y = 0.0;
  }
};

#endif
