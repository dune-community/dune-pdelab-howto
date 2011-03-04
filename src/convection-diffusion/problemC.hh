#ifndef DUNE_PARSOLVE_PROBLEMC_HH
#define DUNE_PARSOLVE_PROBLEMC_HH

#include<math.h>

// function for defining the diffusion tensor
template<typename GV, typename RF>
class k_C
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> >,
      k_C<GV,RF> >
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

public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,k_C<GV,RF> > BaseT;

  k_C (const GV& gv_) : gv(gv_)
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

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  { 
	typename Dune::FieldVector<typename Traits::DomainFieldType,GV::dimension> xglobal;
	xglobal = e.geometry().global(x);

    int ix,iy,iz;

    ix=((int)floor(xglobal[0]/width))%2;
    iy=((int)floor(xglobal[1]/width))%2;
	if (GV::dimension>2)
	  iz=((int)floor(xglobal[2]/width))%2;
	else
	  iz=0;

	RF k;
    if ( iz==0 && iy==0 && ix==0 ) k=K000;
    if ( iz==0 && iy==0 && ix==1 ) k=K001;
    if ( iz==0 && iy==1 && ix==0 ) k=K010;
    if ( iz==0 && iy==1 && ix==1 ) k=K011;
    if ( iz==1 && iy==0 && ix==0 ) k=K100;
    if ( iz==1 && iy==0 && ix==1 ) k=K101;
    if ( iz==1 && iy==1 && ix==0 ) k=K110;
    if ( iz==1 && iy==1 && ix==1 ) k=K111;
	y = k;
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
class K_C
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_C<GV,RF> >
{
  const GV& gv;
  RF K000;
  RF K001;
  RF K010;
  RF K011;
  RF K100;
  RF K101;
  RF K110;
  RF K111;
  RF width;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_C<GV,RF> > BaseT;

  K_C (const GV& gv_) : gv(gv_) 
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

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  { 
	typename Dune::FieldVector<typename Traits::DomainFieldType,GV::dimension> xglobal;
	xglobal = e.geometry().global(x);

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
          y[i][i] = k;
        else
          y[i][j] = 0.0;
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class A0_C
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  A0_C<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,A0_C<GV,RF> > BaseT;

  A0_C (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class F_C
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F_C<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_C<GV,RF> > BaseT;

  F_C (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0; 
  }
};



// constraints parameter class for selecting boundary condition type 
class BCTypeParam_C
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
class B_C
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<
                                                    GV,Dune::PDELab::DiffusionBoundaryCondition::Type,1,
                                                    Dune::FieldVector<
                                                      Dune::PDELab::DiffusionBoundaryCondition::Type,1> >,
                                                  B_C<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::DiffusionBoundaryCondition BC;
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,
     Dune::PDELab::DiffusionBoundaryCondition::Type,1,
     Dune::FieldVector<Dune::PDELab::DiffusionBoundaryCondition::Type,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_C<GV> > BaseT;

  B_C (const GV& gv_) : gv(gv_) {}

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
class G_C
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G_C<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G_C<GV,RF> > BaseT;

  G_C (const GV& gv) : BaseT(gv) {}
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
class J_C
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J_C<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J_C<GV,RF> > BaseT;

  J_C (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0;
	return;
  }
};

// flux as velocity field for the mixed method
template<typename GV, typename RF>
class V_C
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
													  V_C<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_C<GV,RF> > BaseT;

  V_C (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
    y = 0.0;
  }
};


#endif
