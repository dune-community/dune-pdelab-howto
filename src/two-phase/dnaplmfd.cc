// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Solve two-phase flow in porous media with mimetic finite-difference method in pressure-saturation formulation
*/
#include "config.h"
#include<iostream>
#include<vector>
#include<map>
#include <fenv.h>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/superlu.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/finiteelementmap/mimeticfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/intersectionindexset.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/localoperator/twophaseccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/newton/newton.hh>

#include"twophasemfd.hh"

//==============================================================================
// Problem definition
//==============================================================================

const double height = 0.6;
const double width = 1.0;
const double pentry = 755.0;
const double pentry_lens = 1163.0;
const double lens_width_min = 0.3;
const double lens_width_max = 0.7;
const double lens_height_min = 0.2;
const double lens_height_max = 0.3;

// parameter class for local operator
template<typename GV, typename RF>
class TwoPhaseParameter
  : public Dune::PDELab::TwoPhaseParameterInterface<Dune::PDELab::TwoPhaseFullTensorParameterTraits<GV,RF>,
                                                    TwoPhaseParameter<GV,RF> >
{
  static const RF eps1;
  static const RF eps2;

public:
  typedef Dune::PDELab::TwoPhaseFullTensorParameterTraits<GV,RF> Traits;
  enum {dim=GV::Grid::dimension};

  //! constructor
  TwoPhaseParameter ()
  {
    gvector=0; gvector[dim-1]=-9.81;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        if (i == j)
          {
            K_lens[i][j] = 3.32E-11;
            K[i][j] = 6.64E-11;
          }
        else
          {
            K_lens[i][j] = 0.0;
            K[i][j] = 0.0;
          }
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
    return - pentry/pow(1.0-s_l,1/2.5);
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,Traits::GridViewType::Grid::dimension>
      global = e.geometry().global(x);
    for (int i=0; i<dim-1; i++)
      if (global[i]<lens_width_min || global[i]>lens_width_max)
        return pentry/pow(s_l,1/2.5);
    if (global[dim-1]<lens_height_min || global[dim-1]>lens_height_max)
      return pentry/pow(s_l,1/2.5);
    return pentry_lens/sqrt(s_l);
  }

  //! inverse capillary pressure function
  typename Traits::RangeFieldType
  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType pc) const
  {
    return 1.0-(pentry/pc)*(pentry/pc)*sqrt(-pentry/pc);
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,Traits::GridViewType::Grid::dimension>
      global = e.geometry().global(x);
    for (int i=0; i<dim-1; i++)
      if (global[i]<lens_width_min || global[i]>lens_width_max)
        return (pentry/pc)*(pentry/pc)*sqrt(pentry/pc);
    if (global[dim-1]<lens_height_min || global[dim-1]>lens_height_max)
      return (pentry/pc)*(pentry/pc)*sqrt(pentry/pc);
    return (pentry_lens/pc)*(pentry_lens/pc);
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
    return 0.9E-3;
  }

  //! gas phase viscosity
  typename Traits::RangeFieldType
  mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_g) const
  {
    return 1E-3;
  }

  //! absolute permeability
  typename Traits::PermTensorType
  k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return K;
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,Traits::GridViewType::Grid::dimension>
      global = e.geometry().global(x);
    for (int i=0; i<dim-1; i++)
      if (global[i]<lens_width_min || global[i]>lens_width_max)
        return K;
    if (global[dim-1]<lens_height_min || global[dim-1]>lens_height_max)
        return K;
    return K_lens;
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
    return 1460.0;
  }

  //! gas phase molar density
  typename Traits::RangeFieldType
  nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_g) const
  {
    return 1000.0;
  }

  //! liquid phase mass density
  typename Traits::RangeFieldType
  rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType p_l) const
  {
    return 1460.0;
  }

  //! gas phase mass density
  typename Traits::RangeFieldType
  rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType p_g) const
  {
    return 1000.0;
  }

  //! liquid phase boundary condition type
  int
  bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);

    for (int i=0; i<dim-1; i++)
      if (global[i]<eps2 || global[i]>width-eps2)
        return 1; // left & right boundary Dirichlet

    if (global[dim-1]>height-eps2)
      return 0; // top boundary Neumann
    if (global[dim-1]<eps2)
      return 0; // bottom boundary Neumann

    assert(0);
    return -1; // unknown
  }

  //! gas phase boundary condition type
  int
  bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);

    for (int i=0; i<dim-1; i++)
      if (global[i]<eps2 || global[i]>width-eps2)
        return 1; // left & right boundary Dirichlet

    if (global[dim-1]>height-eps2)
      return 0; // top boundary Neumann
    if (global[dim-1]<eps2)
      return 0; // bottom boundary Neumann

    assert(0);
    return -1; // unknown
  }

  //! liquid phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);
    return (height-global[dim-1])*9810.0-pc(*is.inside(),is.geometryInInside().global(x),0.);
  }

  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);
    return (height-global[dim-1])*9810.0;
  }

  //! liquid phase Neumann boundary condition
  typename Traits::RangeFieldType
  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);

    if (global[dim-1]<height-eps2)
      return 0.0;

    for (int i=0; i<dim-1; i++)
      if (global[i]<0.4 || global[i]>0.6)
        return 0.0;

    return -0.075 / 1460.0;
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
  typename Traits::PermTensorType K_lens, K;
};

// Initialize static members. Has to be done out of clas
template<typename GV, typename RF>
const RF TwoPhaseParameter<GV,RF>::eps1 = 1E-6;

template<typename GV, typename RF>
const RF TwoPhaseParameter<GV,RF>::eps2 = 1E-5;

//==============================================================================
// initial conditions for s_w and p_n
//==============================================================================

template<typename GV, typename RF>
class P_n
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> >,
                                          P_n<GV,RF> >
{
  const GV& gv;
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,GV::Grid::dimension,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,P_n<GV,RF> > BaseT;

  P_n (const GV& gv_, const TwoPhaseParameter<GV,RF>& tp_) : gv(gv_), tp(tp_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    Dune::FieldVector<typename Traits::GridViewType::Grid::ctype,dim>
      global = e.geometry().global(x);
    y = (height-global[dim-1])*9810;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }
};

template<typename GV, typename RF>
class S_w
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  S_w<GV,RF> >
{
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,S_w<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  S_w (const GV& gv, const TwoPhaseParameter<GV,RF>& tp_) : BaseT(gv), tp(tp_) {}


  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 0.;
  }
};


//==============================================================================
// Function classes to compute S_n and p_w
//==============================================================================

template<typename  TP, typename PN, typename SW>
class S_n
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PN::Traits::GridViewType,
                                                                           typename PN::Traits::RangeFieldType,PN::Traits::dimRange,
                                                                           typename PN::Traits::RangeType>,
                                          S_n<TP,PN,SW> >
{
  const TP& tp;
  const PN& pn;
  const SW& sw;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PN::Traits::GridViewType,
    typename PN::Traits::RangeFieldType,PN::Traits::dimRange,
    typename PN::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_n<TP,PN,SW> > BaseT;

  S_n (const TP& tp_, const PN& pn_, const SW& sw_) : tp(tp_), pn(pn_), sw(sw_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PN::Traits::RangeType sw_value;
    sw.evaluate(e,x,sw_value);
    y = 1.0 - sw_value;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return pn.getGridView();
  }
};

template<typename  TP, typename PN, typename SW>
class P_w
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PN::Traits::GridViewType,
                                                                           typename PN::Traits::RangeFieldType,PN::Traits::dimRange,
                                                                           typename PN::Traits::RangeType>,
                                          P_w<TP,PN,SW> >
{
  const TP& tp;
  const PN& pn;
  const SW& sw;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PN::Traits::GridViewType,
    typename PN::Traits::RangeFieldType,PN::Traits::dimRange,
    typename PN::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,P_w<TP,PN,SW> > BaseT;

  P_w (const TP& tp_, const PN& pn_, const SW& sw_) : tp(tp_), pn(pn_), sw(sw_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PN::Traits::RangeType pn_value, sw_value;
    pn.evaluate(e,x,pn_value);
    sw.evaluate(e,x,sw_value);
    y = pn_value - tp.pc(e,x,sw_value);
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return pn.getGridView();
  }
};

//==============================================================================
// LinearSolver wrapper class for Newton's method
//==============================================================================

class LinearSolver
{
public:
     explicit LinearSolver(unsigned maxiter_=5000, bool verbose_=true)
         : maxiter(maxiter_), verbose(verbose_)
    {}

    template<class V>
    typename V::ElementType norm(const V& v) const
    {
        return v.two_norm();
    }

    template<class M, class V, class W>
    void apply(M& A, V& z, W& r, typename W::ElementType reduction)
    {
        // Dune::MatrixAdapter<M,V,W> opa(A);
        // Dune::SeqSSOR<M,V,W> ssor(A, 3, 1.0);
        // Dune::BiCGSTABSolver<V> solver(opa, ssor, reduction, maxiter, verbose);
        // Dune::InverseOperatorResult stat;
        // solver.apply(z, r, stat);

        // Dune::printmatrix(std::cout,A.base(),"global stiffness matrix","row",9,1);
        // double max = 0.0;
        // for (int i = 0; i < r.size(); ++i)
        //   {
        //     if (r[i][0] > max)
        //       max = r[i][0];
        //     if (r[i][1] > max)
        //       max = r[i][1];
        //   }

        // for (int i = 0; i < r.size(); ++i)
        //   {
        //     if (std::abs(r[i][0]) > 0.8*max)
        //       std::cout << i << " " << 0 << std::endl;
        //     if (std::abs(r[i][1]) > 0.8*max)
        //       std::cout << i << " " << 1 << std::endl;
        //     // std::cout << std::setw(10) << std::setprecision(2) << std::scientific
        //     //           << r[i][0];
        //     // std::cout << std::setw(10) << std::setprecision(2) << std::scientific
        //     //           << r[i][1];
        //   }

        Dune::SuperLU<typename M::BaseT> superlu(A.base());
        Dune::InverseOperatorResult stat;
        superlu.apply(z, r, reduction, stat);

        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
    }

    const Dune::PDELab::LinearSolverResult<double>& result() const
    {
        return res;
    }

private:
    Dune::PDELab::LinearSolverResult<double> res;
    unsigned maxiter;
    bool verbose;
};

//==============================================================================
// driver
//==============================================================================

template<class GV>
void test (const GV& gv, int timesteps, double timestep)
{
  // some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  int rank = gv.comm().rank();

  // set up index set for intersections
  typedef Dune::PDELab::IntersectionIndexSet<GV> IIS;
  IIS iis(gv);

  // make finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> CellFEM;
  CellFEM cell_fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
  typedef Dune::PDELab::MimeticLocalFiniteElementMap<IIS,DF,RF,dim> FaceFEM;
  FaceFEM face_fem(iis, Dune::GeometryType::cube);

  // make function spaces
  typedef typename Dune::PDELab::ISTLVectorBackend<2> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,CellFEM,
    Dune::PDELab::NoConstraints,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> CellGFS;
  CellGFS cell_gfs(gv, cell_fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,FaceFEM,
    Dune::PDELab::NoConstraints,VBE,
    Dune::PDELab::GridFunctionStaticSize<IIS> > FaceGFS;
  FaceGFS face_gfs(gv, face_fem, iis);
  typedef Dune::PDELab::CompositeGridFunctionSpace
    <Dune::PDELab::GridFunctionSpaceLexicographicMapper,CellGFS,FaceGFS> GFS;
  GFS gfs(cell_gfs, face_gfs);

  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,
    Dune::PDELab::GridFunctionSpaceBlockwiseMapper> TPGFS;
  TPGFS tpgfs(gfs);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // make subspaces (needed for VTK output)
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,0> PhaseN_Subspace;
  PhaseN_Subspace phase_n(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PhaseN_Subspace,0> PhaseN_CellSubspace;
  PhaseN_CellSubspace p_n_gfs(phase_n);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,1> PhaseW_Subspace;
  PhaseW_Subspace phase_w(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<PhaseW_Subspace,0> PhaseW_CellSubspace;
  PhaseW_CellSubspace s_w_gfs(phase_w);

  // make parameter object
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp;

  // initial value function
  typedef P_n<GV,RF> P_nType;
  P_nType p_n_initial(gv,tp);
  typedef S_w<GV,RF> S_wType;
  S_wType s_w_initial(gv,tp);

  // make vector for old time step and initialize
  typedef typename Dune::PDELab::BackendVectorSelector<TPGFS,RF>::Type V;
  V pold(tpgfs);
  pold = 0.0;
  Dune::PDELab::interpolate(p_n_initial,p_n_gfs,pold);
  Dune::PDELab::interpolate(s_w_initial,s_w_gfs,pold);

  // make vector for new time step and initialize
  V pnew(tpgfs);
  pnew = pold;

  // make local operator
  typedef Dune::PDELab::TwoPhaseMFD<TP,V> LOP;
  LOP lop(tp,pold);

  // make grid operator space
  typedef Dune::PDELab::EmptyTransformation C;
  typedef VBE::MatrixBackend MB;
  typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MB,RF,RF,RF,C,C> TPGO;
  TPGO tpgo(tpgfs,tpgfs,lop);

  // initialize face unknowns
  lop.set_time(0.0);
  lop.set_timestep(timestep);
  lop.set_init_mode(true);

  LinearSolver solver(5000, 0);
  Dune::PDELab::Newton<TPGO,LinearSolver,V> newton(tpgo, pnew, solver);
  newton.apply();

  lop.set_init_mode(false);

  // make discrete function objects for pnew and saturations
  typedef Dune::PDELab::DiscreteGridFunction<PhaseN_CellSubspace,V> P_nDGF;
  P_nDGF p_n_dgf(p_n_gfs,pnew);
  typedef Dune::PDELab::DiscreteGridFunction<PhaseW_CellSubspace,V> S_wDGF;
  S_wDGF s_w_dgf(s_w_gfs,pnew);
  typedef S_n<TP,P_nDGF,S_wDGF> S_nDGF;
  S_nDGF s_n_dgf(tp,p_n_dgf,s_w_dgf);
  typedef P_w<TP,P_nDGF,S_wDGF> P_wDGF;
  P_wDGF p_w_dgf(tp,p_n_dgf,s_w_dgf);

  // initialize output of timesteps
  bool graphics = true;
  int filecounter = 0;
  char basename[255];
  sprintf(basename,"dnaplmfd-%01dd",dim);

  // time loop
  RF time = 0.0;
  int k = 1;
  RF total_red = 1.0;
  for (;;)
    {
      if (graphics)
        {
          if (rank==0) std::cout << "writing output file " << filecounter << std::endl;
          Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_nDGF>(p_n_dgf,"p_n"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_wDGF>(p_w_dgf,"p_w"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_nDGF>(s_n_dgf,"s_n"));
          vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_wDGF>(s_w_dgf,"s_w"));
          char fname[255];
          sprintf(fname,"%s-%05d",basename,filecounter);
          vtkwriter.pwrite(fname,"vtk","",Dune::VTK::appendedraw);
          filecounter++;
        }

      if (k > timesteps)
        break;
      lop.set_time(time+timestep);
      lop.set_timestep(timestep);

      // prepare new time step
      if (rank==0) std::cout << "+++ TIME STEP " << k << " tnew=" << time+timestep
                             << " dt=" << timestep << " ++++++++++++++++++" << std::endl;

      LinearSolver solver(5000, 0);
      typedef Dune::PDELab::Newton<TPGO,LinearSolver,V> Newton;
      Newton newton(tpgo, pnew, solver);
      newton.setReassembleThreshold(0.2);
      newton.setAbsoluteLimit(1e-22);
      newton.setReduction(1e-10);
      newton.setMinLinearReduction(1e-9);
      newton.setLineSearchStrategy(Newton::hackbuschReuskenAcceptBest);
      // newton.setMaxIterations(1);
      // newton.setForceIteration(true);
      newton.apply();

      // accept time step
      total_red *= newton.result().reduction;
      // if (total_red < 1e-9)
      //   {
          total_red = 1.0;
          time += timestep;
          pold = pnew;
          ++k;
        // }
    }
}

//==============================================================================
// grid setup
//==============================================================================

int main(int argc, char** argv)
{
  //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
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

    if (argc!=4)
      {
        if(helper.rank()==0)
          std::cout << "usage: ./dnapl <level> <timesteps> <timestep>" << std::endl;
        return 1;
      }

    int maxlevel;
    sscanf(argv[1],"%d",&maxlevel);

    int timesteps;
    sscanf(argv[2],"%d",&timesteps);

    double timestep;
    sscanf(argv[3],"%lg",&timestep);

    // 2D
    if (true)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,2> L; L[0] = width; L[1] = height;
      Dune::FieldVector<int,2> N;    N[0] = 10*(1<<l);   N[1] = 6*(1<<l);
      Dune::FieldVector<bool,2> B(false);
      int overlap=3;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);

      // solve problem :)
      test(grid.leafView(),timesteps,timestep);
    }

    // 3D
    if (false)
    {
      // make grid
      int l=maxlevel;
      Dune::FieldVector<double,3> L; L[0] = width; L[1] = width; L[2] = height;
      Dune::FieldVector<int,3> N;    N[0] = 10*(1<<l);    N[1] = 10*(1<<l);    N[2] = 6*(1<<l);
      Dune::FieldVector<bool,3> B(false);
      int overlap=2;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);

      // solve problem :)
      test(grid.leafView(),timesteps,timestep);
    }

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
