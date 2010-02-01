// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG 
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/diffusion.hh>
#include<dune/pdelab/localoperator/convectiondiffusion.hh>
#include<dune/pdelab/newton/newton.hh>

#include"gridexamples.hh"
#include"problemA.hh"
#include"problemB.hh"
#include"problemC.hh"
#include"problemD.hh"
#include"problemE.hh"
#include"problemF.hh"

//==============================================================================
// Some linear solver variants to be used in Newton's method
//==============================================================================

//! base class for parameter class
template<typename GV, typename RF>
class ConvectionDiffusionProblem : 
  public Dune::PDELab::ConvectionDiffusionParameterInterface<Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>, 
                                                             ConvectionDiffusionProblem<GV,RF> >
{
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! source/reaction term
  typename Traits::RangeFieldType 
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return 0.0;
  }

  //! nonlinearity under gradient
  typename Traits::RangeFieldType 
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return u;
  }

  //! nonlinear scaling of diffusion tensor
  typename Traits::RangeFieldType 
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return 1.0;
  }

  //! tensor permeability
  typename Traits::PermTensorType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType kabs;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        kabs[i][j] = (i==j) ? 1 : 0;
    return kabs;
  }

  //! nonlinear flux vector
  typename Traits::RangeType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    typename Traits::RangeType flux;
    flux[0] = 10 * 1.0 * u*u;
    flux[1] = 10 * 0.5 * u*u;
    return flux;
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
    if (global[0]<1E-6 || global[0]>1-1E-6)
      return 1;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType global = e.geometry().global(x);
    if (global[1]>0.25+global[0]*0.5)
      return 1.0;
    else
      return 0.0;
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType 
  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    return 0.0;
  }
};


//==============================================================================
// Some linear solver variants to be used in Newton's method
//==============================================================================

class SequentialLinearSolver
{
public:
  /*! \brief make a linear solver object
    
    \param[in] maxiter maximum number of iterations to do
    \param[in] verbose print messages if true
  */
  explicit SequentialLinearSolver(unsigned maxiter_=5000, bool verbose_=true)
    : maxiter(maxiter_), verbose(verbose_)
  {}

  /*! \brief compute global norm of a vector
    
    \param[in] v the given vector
  */
  template<class V>
  typename V::ElementType norm(const V& v) const
  {
    return v.two_norm();
  }

  /*! \brief solve the given linear system
    
    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::ElementType reduction)
  {
    Dune::MatrixAdapter<M,V,W> opa(A);
    Dune::SeqSSOR<M,V,W> ssor(A, 3, 1.0);
    Dune::BiCGSTABSolver<V> solver(opa, ssor, reduction, maxiter, verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(z, r, stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
  }
  
  /*! \brief Return access to result data */
  const Dune::PDELab::LinearSolverResult<double>& result() const
  {
    return res;
  }
  
private:
    Dune::PDELab::LinearSolverResult<double> res;
    unsigned maxiter;
    bool verbose;
};


template<class GFS>
class NonoverlappingLinearSolver
{
  typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

public:
  /*! \brief make a linear solver object
    
    \param[in] gfs a grid function space
    \param[in] maxiter maximum number of iterations to do
    \param[in] verbose print messages if true
  */
  explicit NonoverlappingLinearSolver(const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
    : gfs(gfs_), phelper(gfs), maxiter(maxiter_), verbose(verbose_)
  {}

  /*! \brief compute global norm of a vector
    
    \param[in] v the given vector
  */
  template<class V>
  typename V::ElementType norm (const V& v) const
  {
    V x(v); // make a copy because it has to be made consistent
    typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
    PSP psp(gfs,phelper);
    psp.make_consistent(x);
    return psp.norm(x);
  }

  /*! \brief solve the given linear system
    
    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename V::ElementType reduction)
  {
    typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,W> POP;
    POP pop(gfs,A,phelper);
    typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
    PSP psp(gfs,phelper);
    typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,W> PRICH;
    PRICH prich(gfs,phelper);
    int verb=0;
    if (gfs.gridview().comm().rank()==0) verb=verbose;
    Dune::BiCGSTABSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
    Dune::InverseOperatorResult stat;
    solver.apply(z,r,stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
  }

  /*! \brief Return access to result data */
  const Dune::PDELab::LinearSolverResult<double>& result() const
  {
    return res;
  }

private:
  const GFS& gfs;
  PHELPER phelper;
  Dune::PDELab::LinearSolverResult<double> res;
  unsigned maxiter;
  int verbose;
};


template<class GFS, class C>
class OverlappingLinearSolver
{
  typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

public:
  /*! \brief make a linear solver object
    
    \param[in] gfs a grid function space
    \param[in] maxiter maximum number of iterations to do
    \param[in] verbose print messages if true
  */
  explicit OverlappingLinearSolver(const GFS& gfs_, const C& c_, unsigned maxiter_=5000, int verbose_=1)
    : gfs(gfs_), c(c_), phelper(gfs), maxiter(maxiter_), verbose(verbose_)
  {}

  /*! \brief compute global norm of a vector
    
    \param[in] v the given vector
  */
  template<class V>
  typename V::ElementType norm (const V& v) const
  {
    typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> PSP;
    PSP psp(gfs,phelper);
    return psp.norm(v);
  }

  /*! \brief solve the given linear system
    
    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename V::ElementType reduction)
  {
    typedef Dune::PDELab::OverlappingOperator<C,M,V,W> POP;
    POP pop(c,A);
    typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> PSP;
    PSP psp(gfs,phelper);
    typedef Dune::SeqSSOR<M,V,W> SeqPrec;
    SeqPrec seqprec(A,5,1.0);
    typedef Dune::PDELab::OverlappingWrappedPreconditioner<C,GFS,SeqPrec> WPREC;
    WPREC wprec(gfs,seqprec,c,phelper);
    int verb=0;
    if (gfs.gridview().comm().rank()==0) verb=verbose;
    Dune::BiCGSTABSolver<V> solver(pop,psp,wprec,reduction,maxiter,verb);
    Dune::InverseOperatorResult stat;
    solver.apply(z,r,stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
  }

  /*! \brief Return access to result data */
  const Dune::PDELab::LinearSolverResult<double>& result() const
  {
    return res;
  }

private:
  const GFS& gfs;
  const C& c;
  PHELPER phelper;
  Dune::PDELab::LinearSolverResult<double> res;
  unsigned maxiter;
  int verbose;
};


//===============================================================
// solve linear problem
// this is not used; just keep it as a starting point for a
// StationaryLinearSolver class
//===============================================================

template<class GOS, class LS, class V> 
class LinearSolver
{
    typedef typename V::ElementType Real;
    typedef typename GOS::template MatrixContainer<Real>::Type M;
    typedef typename GOS::Traits::TrialGridFunctionSpace::template VectorContainer<Real>::Type W;

public:

  LinearSolver (const GOS& gos_, LS& ls_, V& x_, typename V::ElementType reduction_)
    : gos(gos_), ls(ls_), x(x_), reduction(reduction_)
  {
  }

  void apply ()
  {
    // assemble matrix; optional: assemble only on demand!
    M m(gos); 
    m = 0.0;
    gos.jacobian(x,m);

    // assemble residual
    W r(gos.testGridFunctionSpace(),0.0);
    gos.residual(x,r);  // residual is additive

    // compute correction
    V z(gos.trialGridFunctionSpace(),0.0);
    ls.apply(m,z,r,1E-6); // solver makes right hand side consistent

    // and update
    x -= z;
  }

private:
  const GOS& gos;
  LS& ls;
  V& x;
  typename V::ElementType reduction;
};


template<class GOS, class LS, class V> 
void driver (const GOS& gos, LS& ls, V& x)
{
  Dune::Timer watch;

  // represent operator as a matrix
  typedef typename V::ElementType Real;
  typedef typename GOS::template MatrixContainer<Real>::Type M;
  watch.reset();
  M m(gos); m = 0.0;
  std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
  watch.reset();
  gos.jacobian(x,m); // jacobian is additive except in Dirichlet rows
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;

  // solve the jacobian system
  typedef typename GOS::Traits::TrialGridFunctionSpace::template VectorContainer<Real>::Type W;
  W r(gos.testGridFunctionSpace(),0.0);
  gos.residual(x,r);  // residual is additive
  V z(gos.trialGridFunctionSpace(),0.0);

  // set up parallel solver
  ls.apply(m,z,r,1E-6); // solver makes right hand side consistent
  x -= z;
}

//===============================================================
// Some variants to solve the nonlinear diffusion problem
//===============================================================

// a sequential variant
template<class GV>
void sequential_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_CD<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints(b,gfs,cg);
  std::cout << "constrained dofs=" << cg.size() 
            << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP; 
  LOP lop(param,2);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
  typedef SequentialLinearSolver LS;
  LS ls(5000,1);

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(1);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_sequential_Q1",Dune::VTKOptions::ascii);
}

// a parallel variant for nonoverlapping grids
template<class GV>
void parallel_nonoverlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  int rank = gv.comm().rank();

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem,con);
  con.compute_ghosts(gfs);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  typedef B_A<GV> B;
  B b(gv);
  Dune::PDELab::constraints(b,gfs,cg);
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size() 
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  typedef G_A<GV,Real> G;
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef K_A<GV,Real> K;
  K k(gv);
  typedef A0_A<GV,Real> A0;
  A0 a0(gv);
  typedef F_A<GV,Real> F;
  F f(gv);
  typedef J_A<GV,Real> J;
  J j(gv);
  typedef Dune::PDELab::Diffusion<K,A0,F,B,J> LOP; 
  LOP lop(k,a0,f,b,j,2);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE,true> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
  typedef NonoverlappingLinearSolver<GFS> LS;
  LS ls(gfs,5000,1);

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.apply();
  //driver(gos,ls,x);

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_nonoverlapping_Q1",Dune::VTKOptions::ascii);
}

// a parallel variant for overlapping grids
template<class GV>
void parallel_overlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  int rank = gv.comm().rank();

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  typedef B_A<GV> B;
  B b(gv);
  Dune::PDELab::constraints(b,gfs,cg);
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size() 
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  typedef G_A<GV,Real> G;
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef K_A<GV,Real> K;
  K k(gv);
  typedef A0_A<GV,Real> A0;
  A0 a0(gv);
  typedef F_A<GV,Real> F;
  F f(gv);
  typedef J_A<GV,Real> J;
  J j(gv);
  typedef Dune::PDELab::Diffusion<K,A0,F,B,J> LOP; 
  LOP lop(k,a0,f,b,j,2);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
  typedef OverlappingLinearSolver<GFS,C> LS;
  LS ls(gfs,cg,5000,1);

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.apply();
  //driver(gos,ls,x);

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_overlapping_Q1",Dune::VTKOptions::ascii);
}


//===============================================================
// Main program with grid setup
//===============================================================

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
		  std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }

    // sequential version
    if (helper.size()==1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(8);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(2);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      sequential_Q1(gv);
    }

#if HAVE_MPI
    // nonoverlapping version
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(8);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=0; // needs overlap 0 because overlap elements are not assembled anyway
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      //grid.globalRefine(1);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      parallel_nonoverlapping_Q1(gv);
    }

    // overlapping version
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(8);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=2;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      //grid.globalRefine(1);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      parallel_overlapping_Q1(gv);
    }
#endif

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
