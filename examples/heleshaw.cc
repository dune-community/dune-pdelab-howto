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
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#include"gridexamples.hh"
#include"twophaseop.hh"
#include"permeability_generator.hh"

//==============================================================================
// Problem definition
//==============================================================================

const double height = 0.6;
const double heightw = 0.3; // initial height of water
const double width = 0.4;
const double depth = 0.02;
const double pentry = 1000.0;
const double patm = 1e5;
const double sltop = 0.2;
const double onset = 60.0;
const double period = 60.0;

// parameter class for local operator
template<typename GV, typename RF>
class TwoPhaseParameter 
  : public Dune::PDELab::TwoPhaseParameterInterface<Dune::PDELab::TwoPhaseParameterTraits<GV,RF>, 
                                                    TwoPhaseParameter<GV,RF> >
{
  static const RF eps1 = 1E-6;
  static const RF eps2 = 1E-5;
  
public:
  typedef Dune::PDELab::TwoPhaseParameterTraits<GV,RF> Traits;
  enum {dim=GV::Grid::dimension};

  //! constructor
  TwoPhaseParameter(const GV& gv_) : gv(gv_), is(gv.indexSet()), perm(is.size(0))
  {
    gvector=0; gvector[dim-1]=-9.81;

	typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	typedef typename Traits::DomainFieldType DF;
	double mink=1E100;
	double maxk=-1E100;

    Dune::FieldVector<double,dim> correlation_length(0.4/25.0);
  
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
    return pentry/sqrt(s_l);
  }
	  
  //! inverse capillary pressure function
  typename Traits::RangeFieldType 
  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType pc) const
  {
    typename Traits::RangeFieldType ratio=pentry/pc;
    return ratio*ratio;
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

    // top / bottom
    if (global[dim-1]>height-eps2)
      return 0; // top boundary Neumann
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

    // top / bottom
    if (global[dim-1]>height-eps2)
      return 1; // top boundary Dirichlet
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

    if (global[dim-1]<eps2)
      {
        //        std::cout << "Bingo " << global << std::endl;
        if (time<onset)
          return patm + heightw*9810.0;
        if (fmod(time-onset,period)/period<=0.5)
          return patm + (heightw-0.1)*9810.0;
        else
          return patm + (heightw+0.1)*9810.0;
      }

    return -1; // unknown
  }

  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType 
  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension> 
      global = is.geometry().global(x);

    if (global[dim-1]>height-eps2)
      return patm + pc(*(is.inside()),global,sltop);

    return -1; // unknown
  }

  //! liquid phase Neumann boundary condition
  typename Traits::RangeFieldType 
  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
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

// initial conditions for phase pressures
template<typename GV, typename RF>
class P_l
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_l<GV,RF> >
{
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_l<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  P_l (const GV& gv, const TwoPhaseParameter<GV,RF>& tp_) : BaseT(gv), tp(tp_) {}


  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
      y = patm;
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

      y = patm + tp.pc(e,xlocal,sltop);
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
  const PL pl;
  const PG pg;

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
  const PL pl;
  const PG pg;

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
// driver
//==============================================================================
int rank;

template<class GV> 
void test (const GV& gv, int timesteps, double timestep, double maxtimestep)
{
  // some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType::cube);
  
  // make grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::P0ParallelConstraints,
    Dune::PDELab::ISTLVectorBackend<2>,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,
    Dune::PDELab::GridFunctionSpaceBlockwiseMapper> TPGFS;
  watch.reset();
  Dune::PDELab::P0ParallelConstraints con;
  GFS gfs(gv,fem,con);
  TPGFS tpgfs(gfs);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // make subspaces (needed for VTK output)
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,0> P_lSUB;
  P_lSUB p_lsub(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,1> P_gSUB;
  P_gSUB p_gsub(tpgfs);

  // make parameter object
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp(gv);

  // initial value function
  typedef P_l<GV,RF> P_lType;
  P_lType p_l_initial(gv,tp);
  typedef P_g<GV,RF> P_gType;
  P_gType p_g_initial(gv,tp);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,P_gType> PType;
  PType p_initial(p_l_initial,p_g_initial);

  // make vector for old time step and initialize
  typedef typename TPGFS::template VectorContainer<RF>::Type V;
  V pold(tpgfs);
  Dune::PDELab::interpolate(p_initial,tpgfs,pold);

  // make vector for new time step and initialize
  V pnew(tpgfs);
  pnew = pold;

  // make local operator
  typedef Dune::PDELab::TwoPhaseTwoPointFluxOperator<TP,V> LOP;
  LOP lop(tp,pold);

  // make discrete function objects for pnew and saturations
  typedef Dune::PDELab::DiscreteGridFunction<P_lSUB,V> P_lDGF;
  P_lDGF p_ldgf(p_lsub,pnew);
  typedef Dune::PDELab::DiscreteGridFunction<P_gSUB,V> P_gDGF;
  P_gDGF p_gdgf(p_gsub,pnew);
  typedef S_l<TP,P_lDGF,P_gDGF> S_lDGF; 
  S_lDGF s_ldgf(tp,p_ldgf,p_gdgf);
  typedef S_g<TP,P_lDGF,P_gDGF> S_gDGF; 
  S_gDGF s_gdgf(tp,p_ldgf,p_gdgf);

  // output of timesteps
  bool graphics = true;
  int filecounter = 0;
  char basename[255];
  sprintf(basename,"heleshaw-%01dd",dim);
  if (graphics)
  {
    typedef typename GFS::template VectorContainer<RF>::Type V0;
    V0 partition(gfs,0.0);
    Dune::PDELab::PartitionDataHandle<GFS,V0> pdh(gfs,partition);
    if (gfs.gridview().comm().size()>1)
      gfs.gridview().communicate(pdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V0> DGF0;
    DGF0 pdgf(gfs,partition);

    if (rank==0) std::cout << "writing output file " << filecounter << std::endl;
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_lDGF>(p_ldgf,"p_l"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_gDGF>(p_gdgf,"p_g"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_gDGF>(s_gdgf,"s_g"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(pdgf,"decomposition"));
    char fname[255];
    sprintf(fname,"%s-%05d",basename,filecounter);
    vtkwriter.pwrite(fname,"vtk","",Dune::VTKOptions::binaryappended);
    filecounter++;
  }

  // make constraints map and initialize it from a function
  typedef typename TPGFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(p_initial,tpgfs,cg,false);

  // make grid operator space
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<2,2> MB;
  typedef Dune::PDELab::GridOperatorSpace<TPGFS,TPGFS,LOP,C,C,MB> TPGOS;
  TPGOS tpgos(tpgfs,cg,tpgfs,cg,lop);

  // represent operator as a matrix
  typedef typename TPGOS::template MatrixContainer<RF>::Type M;
  M m(tpgos);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // solver stuff
  typedef Dune::PDELab::ParallelISTLHelper<TPGFS> PHELPER;
  PHELPER phelper(tpgfs);
  typedef Dune::PDELab::OverlappingOperator<C,M,V,V> POP;
  POP pop(cg,m);
  typedef Dune::PDELab::OverlappingScalarProduct<TPGFS,V> PSP;
  PSP psp(tpgfs,phelper);
  int rank = gv.comm().rank();

  // time loop
  RF time = 0.0;
  RF timestepmax=maxtimestep;
  RF timestepscale=1.5;
  for (int k=1; k<=timesteps; k++)
    {
      // prepare new time step
      if (rank==0) std::cout << "+++ TIME STEP " << k << " tnew=" << time+timestep << " dt=" << timestep << std::endl;
      lop.set_time(time+timestep);
      lop.set_timestep(timestep);

      // Newton iteration
      V r(tpgfs,0.0);
      tpgos.residual(pnew,r);
      char buf[64];
      sprintf(buf,"[%02d] ",gv.comm().rank());
      int cols = 14;
      //Dune::printvector(std::cout,r.base(),"initial residual after computation",buf,cols,9,1);
      RF d0 = psp.norm(r);
      RF lastd=d0;
      RF red = 1.0;
      if (rank==0) std::cout << "+++ NEWTON STEP " << 0 << " res=" << d0 << std::endl;
      for (int i=1; i<=25; i++)
        {
          m = 0.0;
          watch.reset();
          tpgos.jacobian(pnew,m);
          if (rank==0) std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;

          typedef Dune::SeqSSOR<M,V,V> SeqPrec;
          SeqPrec seqprec(m,5,1.0);
          typedef Dune::PDELab::OverlappingWrappedPreconditioner<C,TPGFS,SeqPrec> WPREC;
          WPREC  wprec(tpgfs,seqprec,cg,phelper);
//           typedef Dune::PDELab::SuperLUSubdomainSolver<TPGFS,M,V,V> PSUBSOLVE;
//           PSUBSOLVE psubsolve(tpgfs,m);
          int verbose=1;
          if (rank>0) verbose=0;
          Dune::BiCGSTABSolver<V> solver(pop,psp,wprec,
                                         std::max(std::min(1E-3,red*red),1E-10),5000,verbose);
          Dune::InverseOperatorResult stat;  

          V v(tpgfs,0.0);
          solver.apply(v,r,stat);
          //Dune::printvector(std::cout,v.base(),"correction computed in solver",buf,cols,9,1);

          // line search
          bool accept=false;
          RF d;
          for (RF lambda=1.0; lambda>=1E-3; lambda*=0.5)
            {
              V z(tpgfs);
              z = pnew;
              z.axpy(-lambda,v);
              r = 0.0;
              tpgos.residual(z,r);
              d = psp.norm(r);
              if (rank==0) std::cout << "+++ NEWTON line search lambda=" << lambda << " res=" << d << std::endl;
              if (d/lastd<=1-0.25*lambda)
                {
                  pnew = z;
                  red = d/lastd;
                  lastd = d;
                  accept = true;
                  break;
                }
            }

          if (!accept)
            {
              std::cout << "no convergence in line search, exiting" << std::endl;
              exit(1);
            }

          if (rank==0) std::cout << "+++ NEWTON STEP " << i << " res=" << d << std::endl;
          if (d<1E-6*d0) break;
       }

      // accept time step
      time += timestep;
      if (timestep*timestepscale<=timestepmax) timestep*=timestepscale;
      pold = pnew;

      if (graphics)
        {
          if (rank==0) std::cout << "writing output file " << filecounter << std::endl;
          Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_lDGF>(p_ldgf,"p_l"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_gDGF>(p_gdgf,"p_g"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_gDGF>(s_gdgf,"s_g"));
          char fname[255];
          sprintf(fname,"%s-%05d",basename,filecounter);
          vtkwriter.pwrite(fname,"vtk","",Dune::VTKOptions::binaryappended);
          filecounter++;
        }
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
		  std::cout << "usage: ./heleshaw <level> <timesteps> <firsttimestep> <maxtimestep>" << std::endl;
		return 1;
	  }

	int maxlevel;
	sscanf(argv[1],"%d",&maxlevel);

	int timesteps;
	sscanf(argv[2],"%d",&timesteps);

	double timestep;
	sscanf(argv[3],"%lg",&timestep);

	double maxtimestep;
	sscanf(argv[4],"%lg",&maxtimestep);

#if HAVE_MPI
    // 2D
    if (false)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,2> L; L[0] = width; L[1] = height;
      Dune::FieldVector<int,2> N;    N[0] = 40*(1<<l);   N[1] = 60*(1<<l);
      Dune::FieldVector<bool,2> B(false);
      int overlap=3;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
 
      // solve problem :
      test(grid.leafView(),timesteps,timestep,maxtimestep);
    }

    // 3D
    if (true)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,3> L; L[0] = width; L[1] = depth; L[2] = height;
      Dune::FieldVector<int,3> N;    N[0] = 40*(1<<l);    N[1] = 2*(1<<l);    N[2] = 60*(1<<l);
      Dune::FieldVector<bool,3> B(false);
      int overlap=2;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
      
      // solve problem :)
      test(grid.leafView(),timesteps,timestep,maxtimestep);
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
