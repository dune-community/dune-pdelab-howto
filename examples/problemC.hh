#ifndef DUNE_PARSOLVE_PROBLEMC_HH
#define DUNE_PARSOLVE_PROBLEMC_HH

#include<math.h>

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

	RF k;
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
  
  inline void fillPermArray(std::vector<RF>& perm) const
  {
    // This is not really congruent to evaluate as there
    // the permeability at the dixretization point is taken
    // in contrast here the value in the element center is taken.
    perm.resize(gv.indexSet().size(0));
    
    typename Traits::RangeType y;
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
	int id = gv.indexSet().index(*it);
        Dune::GeometryType gt = it->geometry().type();
	typedef typename Traits::DomainFieldType DF;
	const int dim = GV::dimension;
        typename Traits::DomainType localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
	evaluate(*it, localcenter,y);
	perm[id]=y[0][0];
      }
    
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
	y = 1; 
  }
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B_C
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<GV,int,1,
                                                                             Dune::FieldVector<int,1> >,
                                                  B_C<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_C<GV> > BaseT;

  B_C (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    y = 1; // Dirichlet
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
	y = 0.0; 
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
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_C<GV,RF> > BaseT;

  V_C (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
    y = 0.0;
  }
};


#endif
