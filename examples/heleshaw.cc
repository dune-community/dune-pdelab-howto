// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p0constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/localoperator/twophaseccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include"gridexamples.hh"
#include"componenttransportop.hh"

#undef RANDOMPERM
#ifdef RANDOMPERM
#include"permeability_generator.hh"
#endif

//==============================================================================
// Problem definition
//==============================================================================

const double height = 0.6;
const double heightw = 0.2; // initial height of water
const double width = 0.4;
const double depth = 0.02;
const double pentry = 1000.0;
const double patm = 1e5;
const double onset = 6000000.0;
const double period = 60.0;

// parameter class for local operator
template<typename GV, typename RF>
class TwoPhaseParameter 
  : public Dune::PDELab::TwoPhaseParameterInterface<Dune::PDELab::TwoPhaseParameterTraits<GV,RF>, 
                                                    TwoPhaseParameter<GV,RF> >
{
  static const RF eps1 = 1E-4; // regularization in kr
  static const RF eps2 = 1E-5; // eps in geometry calculations
  
public:
  typedef Dune::PDELab::TwoPhaseParameterTraits<GV,RF> Traits;
  enum {dim=GV::Grid::dimension};

  //! constructor
  TwoPhaseParameter(const GV& gv_) : gv(gv_), is(gv.indexSet()), perm(is.size(0))
  {
    gvector=0; gvector[dim-1]=-9.81;

	typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	typedef typename Traits::DomainFieldType DF;
  
#ifdef RANDOMPERM
	double mink=1E100;
	double maxk=-1E100;
    Dune::FieldVector<double,dim> correlation_length(0.4/50.0);
	EberhardPermeabilityGenerator<GV::dimension> field(correlation_length,1,0.0,5000,-1083);
	for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
	  {
		int id = is.index(*it);
        Dune::GeometryType gt = it->geometry().type();
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
		perm[id]=field.eval(globalcenter);
		mink = std::min(mink,perm[id]);
		maxk = std::max(maxk,perm[id]);
	  }
	std::cout << "       mink=" << mink << "               maxk=" << maxk << std::endl;
	std::cout << "log10(mink)=" << log10(mink) << " log10(maxk)=" << log10(maxk) << std::endl;
#else
	for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
	  {
		int id = is.index(*it);
		perm[id]=1.0;
	  }
#endif
  }

  //! porosity
  typename Traits::RangeFieldType 
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.4;
  }

  //! capillary pressure function
  typename Traits::RangeFieldType 
  pc (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType s_l) const
  {
    const typename Traits::RangeFieldType epsreg=1e-5;
    const typename Traits::RangeFieldType S0reg=1.0-1e-5;
    if (s_l>S0reg-epsreg)
      return (S0reg-s_l)*pentry/epsreg;
    else
      return pentry/sqrt(s_l/(S0reg-epsreg));
  }
	  
  //! inverse capillary pressure function
  typename Traits::RangeFieldType 
  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType pc) const
  {
    typename Traits::RangeFieldType ratio=pentry/pc;
    const typename Traits::RangeFieldType epsreg=1e-5;
    const typename Traits::RangeFieldType S0reg=1.0-1e-5;
    if (pc<=pentry)
      return S0reg-epsreg/ratio;
    else
      return (S0reg-epsreg)*ratio*ratio;
  }
	  
  //! liquid phase relative permeability
  typename Traits::RangeFieldType 
  kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType s_l) const
  {
    if (s_l<=eps1) return 0.0; else return (s_l-eps1)*(s_l-eps1);
  }

  //! gas phase relative permeability
  typename Traits::RangeFieldType 
  kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType s_g) const
  {
    if (s_g<=eps1) return 0.0; else return (s_g-eps1)*(s_g-eps1);
  }

  //! liquid phase viscosity
  typename Traits::RangeFieldType 
  mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_l) const
  {
    return 1E-3;
  }

  //! gas phase viscosity
  typename Traits::RangeFieldType 
  mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_g) const
  {
    return 4.65E-5;
  }
	  
  //! absolute permeability (scalar!)
  typename Traits::RangeFieldType 
  k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    int id = is.index(e);
    return perm[id]*6.64E-11;
  }

  //! gravity vector
  const typename Traits::RangeType& gravity () const
  {
    return gvector;
  }

  //! liquid phase molar density
  typename Traits::RangeFieldType 
  nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_l) const
  {
    return 1000.0;
  }

  //! gas phase molar density
  typename Traits::RangeFieldType 
  nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_g) const
  {
    return p_g/(287.2*300.0);
  }

  //! liquid phase mass density
  typename Traits::RangeFieldType 
  rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
         typename Traits::RangeFieldType p_l) const
  {
    return 1000.0;
  }

  //! gas phase mass density
  typename Traits::RangeFieldType 
  rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
         typename Traits::RangeFieldType p_g) const
  {
    return p_g/(287.2*300.0);
  }
	  
  //! liquid phase boundary condition type
  int
  bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    // left / right
    if (global[0]<eps2 || global[0]>width-eps2)
      return 0; // left & right boundary: Neumann

    // top 
    if (global[dim-1]>height-eps2)
      return 0; // top boundary Neumann

    // bottom
    if (global[dim-1]<eps2)
      return 1; // bottom boundary Dirichlet

    // front / back
    if (dim==3)
      {
        if (global[1]<eps2 || global[1]>depth-eps2)
          return 0; // left & right boundary: Neumann
     }

    return -1; // unknown
  }

  //! gas phase boundary condition type
  int
  bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    // left / right
    if (global[0]<eps2 || global[0]>width-eps2)
      return 0; // left & right boundary: Neumann

    // top
    if (global[dim-1]>height-eps2)
      {
        // { 17.II.2010: influx at top
        {
          typename Traits::RangeFieldType w = 0.03;        
          if (global[0]>0.1-w && global[0]<0.1+w)
            return 0;
          if (global[0]>0.2-w && global[0]<0.2+w)
            return 0;
          if (global[0]>0.3-w && global[0]<0.3+w)
            return 0;
        }
        // }

        return 1; // top boundary Dirichlet
      }

    // bottom
    if (global[dim-1]<eps2)
      return 0; // bottom boundary Neumann

    // front / back
    if (dim==3)
      {
        if (global[1]<eps2 || global[1]>depth-eps2)
          return 0; // left & right boundary: Neumann
     }

    return -1; // unknown
  }

  //! liquid phase Dirichlet boundary condition
  typename Traits::RangeFieldType 
  g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    return  patm - pentry + (heightw-global[dim-1])*9810.0;
  }

  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType 
  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    return patm + (heightw-global[dim-1])*9.81*patm/(287.2*300.0);
  }

  //! liquid phase Neumann boundary condition
  typename Traits::RangeFieldType 
  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    // { 17.II.2010: influx at top
    if (global[dim-1]>height-eps2)
      {
        typename Traits::RangeFieldType flux = -0.075;
        if (time<1e-6) 
          flux = 0.0;
        else
          {
            if (time>0.0)
              flux = -0.075;
            else
              flux = -0.075*(time/20.0);
          }
        typename Traits::RangeFieldType w = 0.03;        
        if (global[0]>0.1-w && global[0]<0.1+w)
          return flux;
        if (global[0]>0.2-w && global[0]<0.2+w)
          return flux;
        if (global[0]>0.3-w && global[0]<0.3+w)
          return flux;
      }
    // }

    return 0.0;
  }

  //! gas phase Neumann boundary condition
  typename Traits::RangeFieldType 
  j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  //! liquid phase source term
  typename Traits::RangeFieldType 
  q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }
  
  //! gas phase source term
  typename Traits::RangeFieldType 
  q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

private:
  typename Traits::RangeType gvector;
  const GV& gv;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;
};


template<typename GV, typename RF>
class P_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> >,
                                          P_l<GV,RF> >
{
  const GV& gv;
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,P_l<GV,RF> > BaseT;

  P_l (const GV& gv_, const TwoPhaseParameter<GV,RF>& tp_) : gv(gv_), tp(tp_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,dim> 
      x = e.geometry().global(xlocal);

    y = patm - pentry + (heightw-x[dim-1])*9810.0;
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }
};


template<typename GV, typename RF>
class P_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> >,
                                          P_g<GV,RF> >
{
  const GV& gv;
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,P_g<GV,RF> > BaseT;

  P_g (const GV& gv_, const TwoPhaseParameter<GV,RF>& tp_) : gv(gv_), tp(tp_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,dim> 
      x = e.geometry().global(xlocal);

    y =  patm + (heightw-x[dim-1])*9.81*patm/(287.2*300.0);
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }
};



//==============================================================================
// saturation output
//==============================================================================

template<typename  T, typename PL, typename PG>
class S_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_l<T,PL,PG> >
{
  const T& t;
  const PL& pl;
  const PG& pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
    typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
    typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_l<T,PL,PG> > BaseT;

  S_l (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = t.s_l(e,x,pg_value-pl_value);
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }
};

template<typename  T, typename PL, typename PG>
class S_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_g<T,PL,PG> >
{
  const T& t;
  const PL& pl;
  const PG& pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
    typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
    typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_g<T,PL,PG> > BaseT;

  S_g (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = 1.0-t.s_l(e,x,pg_value-pl_value);
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }
};


//==============================================================================
// Problem definition : component transport
//==============================================================================


//! base class for parameter class
template<typename GV, typename RF, typename SOLDDGF, typename SNEWDGF, typename UDGF>
class Transport1 : 
  public Dune::PDELab::ModifiedTransportSpatialParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>, 
                                                                  Transport1<GV,RF,SOLDDGF,SNEWDGF,UDGF> >,
  public Dune::PDELab::ModifiedTransportTemporalParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>, 
                                                                  Transport1<GV,RF,SOLDDGF,SNEWDGF,UDGF> >
{
  typedef TwoPhaseParameter<GV,RF> TP;
  enum {dim=GV::Grid::dimension};

public:
  typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;

  Transport1 (const TP& tp_, const SOLDDGF& solddgf_, const SNEWDGF& snewdgf_, const UDGF& udgf_)
    : tp(tp_), solddgf(solddgf_), snewdgf(snewdgf_), udgf(udgf_)
  {}

  //! capacity function
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType s_new;
    typename SOLDDGF::Traits::RangeType s_old;
    snewdgf.evaluate(e,x,s_new);
    solddgf.evaluate(e,x,s_old);
    RF factor = (time-t0)/dt;
    return tp.nu_l(e,x,0.0)*tp.phi(e,x)*( (1.0-factor)*s_old + factor*s_new );
  }

  //! saturation at new time level
  typename Traits::RangeFieldType 
  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType y;
    snewdgf.evaluate(e,x,y);
    return y;
  }

  //! saturation at old time level
  typename Traits::RangeFieldType 
  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SOLDDGF::Traits::RangeType y;
    solddgf.evaluate(e,x,y);
    return y;
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType velo;
    udgf.evaluate(e,x,velo);
    return velo;
  }

  //! tensor permeability
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source/reaction term
  typename Traits::RangeFieldType 
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 0.0;
  }

  //! boundary condition type function
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  int
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::RangeType global = is.geometry().global(x);
    return 1;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType global = e.geometry().global(x);

    if (time<500 || time>600) return 0;
    if (global[dim-1]<height-1E-7) return 0;

    if (global[0]>=0.19 && global[0]<=0.21)
        return 1.0;

    return 0.0;
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType 
  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

  //! to be called once before each time step
  void preStep (RF time_, RF dt_, int stages)
  {
    t0 = time_;
    dt = dt_;
  }
      
private:
  const TP& tp; // store access to two phase parameters
  const SOLDDGF& solddgf;
  const SNEWDGF& snewdgf;
  const UDGF& udgf;
  RF time, t0, dt;
};


// initial conditions for component concentration
template<typename GV, typename RF>
class Conc
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  Conc<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Conc<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  Conc (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
 	y = 0;
  }
};


//==============================================================================
// driver
//==============================================================================
int rank;

template<class GV> 
void test (const GV& gv, double Tend, double timestep, double maxtimestep, int level)
{
  // <<<1>>> choose some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType::cube);
  typedef Dune::PDELab::P0ParallelConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<2> VBE2;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE2,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,
    Dune::PDELab::GridFunctionSpaceBlockwiseMapper> TPGFS;
  watch.reset();
  CON con;
  GFS gfs(gv,fem,con);
  TPGFS tpgfs(gfs);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // <<<2b>>> make subspaces for visualization
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,0> P_lSUB;
  P_lSUB p_lsub(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,1> P_gSUB;
  P_gSUB p_gsub(tpgfs);

  // make parameter object
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp(gv);

  // <<<4>>> initial value function
  typedef P_l<GV,RF> P_lType;
  P_lType p_l_initial(gv,tp);
  typedef P_g<GV,RF> P_gType;
  P_gType p_g_initial(gv,tp);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,P_gType> PType;
  PType p_initial(p_l_initial,p_g_initial);

  // <<<5>>> make vector for old time step and initialize
  typedef typename TPGFS::template VectorContainer<RF>::Type V;
  V pold(tpgfs);
  Dune::PDELab::interpolate(p_initial,tpgfs,pold);

  // <<<6>>> make vector for new time step and initialize
  V pnew(tpgfs);
  pnew = pold;

  // <<<7>>> make discrete function objects for pnew and saturations
  typedef Dune::PDELab::DiscreteGridFunction<P_lSUB,V> P_lDGF;
  P_lDGF p_ldgf(p_lsub,pnew);
  typedef Dune::PDELab::DiscreteGridFunction<P_gSUB,V> P_gDGF;
  P_gDGF p_gdgf(p_gsub,pnew);
  typedef S_l<TP,P_lDGF,P_gDGF> S_lDGF; 
  S_lDGF s_ldgf(tp,p_ldgf,p_gdgf);
  typedef S_g<TP,P_lDGF,P_gDGF> S_gDGF; 
  S_gDGF s_gdgf(tp,p_ldgf,p_gdgf);

  P_lDGF pold_ldgf(p_lsub,pold);
  P_gDGF pold_gdgf(p_gsub,pold);
  S_lDGF sold_ldgf(tp,pold_ldgf,pold_gdgf);
  S_gDGF sold_gdgf(tp,pold_ldgf,pold_gdgf);

  // <<<9>>> velocity grid functions
  typedef Dune::PDELab::V_l<TP,P_lDGF,P_gDGF> VliquidDGF;
  VliquidDGF vliquiddgf(tp,p_ldgf,p_gdgf);
  vliquiddgf.set_time(0);
  VliquidDGF vliquidolddgf(tp,pold_ldgf,pold_gdgf);
  vliquidolddgf.set_time(0);
  typedef Dune::PDELab::V_g<TP,P_lDGF,P_gDGF> VgasDGF;
  VgasDGF vgasdgf(tp,p_ldgf,p_gdgf);
  vgasdgf.set_time(0);
  VgasDGF vgasolddgf(tp,pold_ldgf,pold_gdgf);
  vgasolddgf.set_time(0);

  // <<<8>>> make constraints map and initialize it from a function
  typedef typename TPGFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(p_initial,tpgfs,cg,false);

  // <<<9>>> make grid operator space
  typedef Dune::PDELab::TwoPhaseTwoPointFluxOperator<TP> LOP;
  LOP lop(tp);
  typedef Dune::PDELab::TwoPhaseOnePointTemporalOperator<TP> MLOP;
  MLOP mlop(tp);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<2,2> MBE;
  Dune::PDELab::ImplicitEulerParameter<RF> method;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<RF,V,TPGFS,TPGFS,LOP,MLOP,C,C,MBE> IGOS;
  IGOS igos(method,tpgfs,cg,tpgfs,cg,lop,mlop);

  // <<<10>>> Make a linear solver 
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<TPGFS,C> LS;
  LS ls(tpgfs,cg,5000,5,0);

  // <<<11>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<IGOS,LS,V> PDESOLVER;
  PDESOLVER tnewton(igos,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(3);
  tnewton.setReduction(1e-6);
  tnewton.setMinLinearReduction(1e-3);
  tnewton.setAbsoluteLimit(1e-7);

  // <<<12>>> time-stepper
  Dune::PDELab::OneStepMethod<RF,IGOS,PDESOLVER,V,V> osm(method,igos,tnewton);
  osm.setVerbosityLevel(2);

  //========================
  // transport related stuff
  //========================

  // grid function space
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE1;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    VBE1,Dune::PDELab::SimpleGridFunctionStaticSize> CGFS;
  CGFS cgfs(gv,fem); 

  // initial value function
  typedef Conc<GV,RF> ConcType;
  ConcType conc_initial(gv);

  // vectors and initialization
  typedef typename CGFS::template VectorContainer<RF>::Type CV;
  CV cold(cgfs); 
  Dune::PDELab::interpolate(conc_initial,cgfs,cold);
  CV cnew(cgfs);
  cnew = cold;
  
  // concentration grid functions
  typedef Dune::PDELab::DiscreteGridFunction<CGFS,CV> ConcDGF;
  ConcDGF concdgf(cgfs,cnew);

  // concentration constraints
  typedef typename CGFS::template ConstraintsContainer<RF>::Type CC;
  CC ccg;
  ccg.clear();
  Dune::PDELab::constraints(conc_initial,cgfs,ccg,false);

  //  make grid operator space for transport problem
  typedef Transport1<GV,RF,S_lDGF,S_lDGF,VliquidDGF> CTP;
  CTP ctp(tp,sold_ldgf,s_ldgf,vliquiddgf);
  typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP> CLOP; 
  CLOP clop(ctp);
  typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> CSLOP; 
  CSLOP cslop(ctp);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE1;
  Dune::PDELab::ExplicitEulerParameter<RF> cmethod;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<RF,CV,CGFS,CGFS,CLOP,CSLOP,CC,CC,MBE1> CIGOS;
  CIGOS cigos(cmethod,cgfs,ccg,cgfs,ccg,clop,cslop);

  // concentration linear solver
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<CGFS> CLS;
  CLS cls(cgfs);

  // <<<13>>> graphics for initial value
  bool graphics = true;
  char basename[255];
  sprintf(basename,"heleshaw-infiltration-IE-%01dl%01dd",level,dim);
  Dune::PDELab::FilenameHelper fn(basename);
  if (graphics)
  {
//     typedef typename GFS::template VectorContainer<RF>::Type V0;
//     V0 partition(gfs,0.0);
//     Dune::PDELab::PartitionDataHandle<GFS,V0> pdh(gfs,partition);
//     if (gfs.gridview().comm().size()>1)
//       gfs.gridview().communicate(pdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
//     typedef Dune::PDELab::DiscreteGridFunction<GFS,V0> DGF0;
//     DGF0 pdgf(gfs,partition);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_lDGF>(p_ldgf,"p_l"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_gDGF>(p_gdgf,"p_g"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_gDGF>(s_gdgf,"s_g"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VliquidDGF>(vliquiddgf,"liquid velocity"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VgasDGF>(vgasdgf,"gas velocity"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"mass fraction"));
    //    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(pdgf,"decomposition"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<14>>> time loop
  RF time = 0.0;
  RF timestepmax=maxtimestep;
  RF timestepscale=1.1;
  RF dtmin = 1e-6;
  bool implicit=true;
  while (time<Tend)
    {
      // do time step
      if (time<500.0) // assume it is stationary ...
        {
          try {
            osm.apply(time,timestep,pold,pnew);
          }
          catch (Dune::PDELab::NewtonLineSearchError) {
            timestep *= 0.5;
            if (timestep<dtmin) throw;
            continue;
          }
        }
      else timestep = timestepmax;

      vliquiddgf.set_time(time+timestep);
      vgasdgf.set_time(time+timestep);

      if (time>=500.0)
        {
          if (implicit)
            {
              // linear solver for transport problem 
              typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<CGFS,CC> CLS;
              CLS cls(cgfs,ccg,5000,5,1);
              
              // make Newton for transport problem
              typedef Dune::PDELab::Newton<CIGOS,CLS,CV,CV> CPDESOLVER;
              CPDESOLVER ctnewton(cigos,cls);
              ctnewton.setReassembleThreshold(0.0);
              ctnewton.setVerbosityLevel(2);
              ctnewton.setReduction(0.9);
              ctnewton.setMinLinearReduction(1e-10);
              
              // <<<16>>> time-steppers
              Dune::PDELab::ImplicitEulerParameter<RF> cmethod;
              cigos.setMethod(cmethod);
              Dune::PDELab::OneStepMethod<RF,CIGOS,CPDESOLVER,CV,CV> cosm(cmethod,cigos,ctnewton);
              cosm.setVerbosityLevel(1);
              
              // transport time stepping
              RF ctime = time;
              RF ctimestep = timestep/1.0;
              while (ctime<time+timestep-1e-8)
                {
                  cosm.apply(ctime,ctimestep,cold,cnew);
                  ctime += ctimestep;
                  cold = cnew;
                }
            }
          else
            {
              // concentration time stepper
              typedef Dune::PDELab::CFLTimeController<RF,CIGOS> TC;
              TC tc(0.999,time+timestep,cigos);
              Dune::PDELab::ExplicitOneStepMethod<RF,CIGOS,CLS,CV,CV,TC> cosm(cmethod,cigos,cls,tc);
              cosm.setVerbosityLevel(2);
              
              // transport time stepping
              RF ctime = time;
              RF ctimestep = timestep;
              while (ctime<time+timestep-1e-8)
                {
                  ctimestep = cosm.apply(ctime,ctimestep,cold,cnew);
                  ctime += ctimestep;
                  cold = cnew;
                }
            }
        }

      // graphical output
      if (graphics)
        {
          Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_lDGF>(p_ldgf,"p_l"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_gDGF>(p_gdgf,"p_g"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_gDGF>(s_gdgf,"s_g"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VliquidDGF>(vliquiddgf,"liquid velocity"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VgasDGF>(vgasdgf,"gas velocity"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"mass fraction"));
          vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
          fn.increment();
        }

      // accept time step
      vliquidolddgf.set_time(time);
      vgasolddgf.set_time(time);
      time += timestep;
      if (timestep*timestepscale<=timestepmax) 
        timestep*=timestepscale;
      else
        timestep = timestepmax;
      pold = pnew;
    }
}

//==============================================================================
// grid setup
//==============================================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
		if(helper.rank()==0)
		  std::cout << "parallel run on " << helper.size() << " processes" << std::endl;
	  }
    rank = helper.rank();

	if (argc!=5)
	  {
		if(helper.rank()==0)
		  std::cout << "usage: ./heleshaw <level> <end time> <firsttimestep> <maxtimestep>" << std::endl;
		return 1;
	  }

	int maxlevel;
	sscanf(argv[1],"%d",&maxlevel);

	double Tend;
	sscanf(argv[2],"%lg",&Tend);

	double timestep;
	sscanf(argv[3],"%lg",&timestep);

	double maxtimestep;
	sscanf(argv[4],"%lg",&maxtimestep);

#if HAVE_MPI
    // 2D
    if (true)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,2> L; L[0] = width; L[1] = height;
      Dune::FieldVector<int,2> N;    N[0] = 40*(1<<l);   N[1] = 60*(1<<l);
      Dune::FieldVector<bool,2> B(false);
      int overlap=3;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
 
      // solve problem :
      test(grid.leafView(),Tend,timestep,maxtimestep,l);
    }

    // 3D
    if (false)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,3> L; L[0] = width; L[1] = depth; L[2] = height;
      Dune::FieldVector<int,3> N;    N[0] = 40*(1<<l);    N[1] = 2*(1<<l);    N[2] = 60*(1<<l);
      Dune::FieldVector<bool,3> B(false);
      int overlap=2;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
      
      // solve problem :)
      test(grid.leafView(),Tend,timestep,maxtimestep,l);
    }
#endif

	// test passed
	return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
} 
