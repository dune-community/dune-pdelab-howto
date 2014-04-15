#ifndef DUNE_PARSOLVE_PROBLEMD_HH
#define DUNE_PARSOLVE_PROBLEMD_HH

#include<math.h>
#include"../utility/permeability_generator.hh"

// function for defining the diffusion tensor
template<typename GV, typename RF>
class k_D
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
                                                                           1,Dune::FieldVector<RF,1> >,
      k_D<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,k_D<GV,RF> > BaseT;

  k_D (const GV& gv_, Dune::FieldVector<double,GV::dimension> correlation_length,
     double variance = 1.0, double mean = 0.0, long modes = 1000, long seed = -1083)
  : gv(gv_), is(gv.indexSet()), perm(is.size(0))
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
          Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
    perm[id]=field.eval(globalcenter);
    mink = std::min(mink,log10(perm[id]));
    maxk = std::max(maxk,log10(perm[id]));
    }
  std::cout << "log10(mink)=" << mink << " log10(maxk)=" << maxk << std::endl;
  }

  k_D ( const GV& gv_, const std::vector<RF>& perm_)
    : gv(gv_), is(gv.indexSet()), perm(perm_)
  {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
  y = perm[is.index(e)];
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

private:
  const GV& gv;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;
};

// function for defining the diffusion tensor
template<typename GV, typename RF>
class K_D
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_D<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_D<GV,RF> > BaseT;

  K_D (const GV& gv_, Dune::FieldVector<double,GV::dimension> correlation_length,
     double variance = 1.0, double mean = 0.0, long modes = 1000, long seed = -1083)
  : gv(gv_), is(gv.indexSet()), perm(is.size(0))
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
          Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
    perm[id]=field.eval(globalcenter);
    perm[id]=1.0;
    mink = std::min(mink,log10(perm[id]));
    maxk = std::max(maxk,log10(perm[id]));
    }
  std::cout << "log10(mink)=" << mink << " log10(maxk)=" << maxk << std::endl;
  }

  K_D ( const GV& gv_, const std::vector<RF>& perm_)
    : gv(gv_), is(gv.indexSet()), perm(perm_)
  {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
      y[i][i] = perm[is.index(e)];
        else
          y[i][j] = 0.0;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

  inline const RF& getElementPermeability(const typename GV::template Codim<0>::EntityPointer& e) const
  {
    return perm[is.index(*e)];
  }

private:
  const GV& gv;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;
};

// function for defining the source term
template<typename GV, typename RF>
class A0_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  A0_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,A0_D<GV,RF> > BaseT;

  A0_D (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class F_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_D<GV,RF> > BaseT;

  F_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
  y = 0;
  }
};


// constraints parameter class for selecting boundary condition type
class BCTypeParam_D
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
    return false;  // Dirichlet
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
class B_D
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<
                                                    GV,Dune::PDELab::DiffusionBoundaryCondition::Type,1,
                                                    Dune::FieldVector<
                                                      Dune::PDELab::DiffusionBoundaryCondition::Type,1> >,
                                                  B_D<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::DiffusionBoundaryCondition BC;
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,
     Dune::PDELab::DiffusionBoundaryCondition::Type,1,
     Dune::FieldVector<Dune::PDELab::DiffusionBoundaryCondition::Type,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_D<GV> > BaseT;

  B_D (const GV& gv_) : gv(gv_) {}

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
class G_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G_D<GV,RF> > BaseT;

  G_D (const GV& gv) : BaseT(gv) {}
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
class J_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J_D<GV,RF> > BaseT;

  J_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
  y = 0;
  return;
  }
};

// flux as velocity field for the mixed method
template<typename GV, typename RF>
class V_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                            V_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_D<GV,RF> > BaseT;

  V_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};


#endif
