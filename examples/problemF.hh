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

// boundary grid function selecting boundary conditions 
template<typename GV>
class B_F
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<GV,int,1,
                                                                             Dune::FieldVector<int,1> >,
                                                  B_F<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
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
	  y = 1; // Dirichlet
	else
	  y = 0;
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

#endif
