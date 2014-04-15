// -*- tab-width: 4; indent-tabs-mode: nil -*-
// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve heat equation with high order in space and time (known analytic solution)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
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
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/convectiondiffusion.hh>
#include<dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include"../utility/gridexamples.hh"
#include"l2interpolationerror.hh"

//==============================================================================
// Parameter class for the convection diffusion problem
//==============================================================================

const double pi = 3.141592653589793238462643;

// grid function for analytic solution at T=0.125
template<typename GV, typename RF>
class U : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  U<GV,RF> > {
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,U<GV,RF> > B;

  U (const GV& gv) : B(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = sin(2.0*pi*0.125) * sin(3.0*pi*x[0]) * sin(2.0*pi*x[1]);
  }
};

//! base class for parameter class
template<typename GV, typename RF>
class ConvectionDiffusionProblem :
  public Dune::PDELab::ConvectionDiffusionParameterInterface<
  Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>,
  ConvectionDiffusionProblem<GV,RF>
  >
{
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! source/reaction term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    typename Traits::RangeType global = e.geometry().global(x);
    typename Traits::RangeFieldType X = sin(3.0*pi*global[0]);
    typename Traits::RangeFieldType Y = sin(2.0*pi*global[1]);
    return X*Y*(2.0*pi*cos(2.0*pi*time)+13.0*pi*pi*sin(2.0*pi*time));
    // exact solution is u(x,y,t) = sin(2*pi*t) * sin(3*pi*x) * sin(2*pi*y)
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
    flux[0] = 0.0;
    flux[1] = 0.0;
    return flux;
  }

  //! boundary condition type function
  template<typename I>
  bool isDirichlet(
                   const I & intersection,               /*@\label{bcp:name}@*/
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {

    //Dune::FieldVector<typename I::ctype, I::dimension>
    //  xg = intersection.geometry().global( coord );

    return true;  // Dirichlet b.c. on all boundaries
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
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

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  RF time;
};

//===============================================================
// Some variants to solve the nonlinear diffusion problem
//===============================================================

// a sequential variant
template<class GV>
void sequential (const GV& gv, int t_level)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;

  // <<<2>>> Make grid function space
  const int degree=2;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::BCTypeParam_CD<Param> bctype(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints( bctype, gfs, cg );
  std::cout << "constrained dofs=" << cg.size()
            << " of " << gfs.globalSize() << std::endl;

  // <<<5>>> Make grid operator space for time-dependent problem
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP;
  LOP lop(param,4);
  typedef Dune::PDELab::L2 MLOP;
  MLOP mlop(4);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  //Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::Alexander3Parameter<Real> method;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO0;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,Real,Real,Real,C,C> GO1;
  GO1 go1(gfs,cg,gfs,cg,mlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  typedef typename IGO::Traits::Domain V;

  // <<<6>>> Make a linear solver
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,0);

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
  PDESOLVER tnewton(igo,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(0);
  tnewton.setReduction(0.9);
  tnewton.setMinLinearReduction(1e-9);

  // <<<8>>> time-stepper
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,V,V> osm(method,igo,tnewton);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("instationarytest_Q1");
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
    fn.increment();
  }

  // <<<11>>> time loop
  Real time = 0.0;
  int N=1; for (int i=0; i<t_level; i++) N *= 2;
  Real dt = 0.125/N;
  V x(gfs,0.0);
  param.setTime(dt);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);
  for (int i=1; i<=N; i++)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,x);
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
      vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
      fn.increment();

      // advance time step
      //       std::cout.precision(8);
      //       std::cout << "solution maximum: "
      //                 << std::scientific << x.infinity_norm() << std::endl;
      xold = x;
      time += dt;
    }

  // evaluate discretization error
  U<GV,Real> u(gv);
  std::cout.precision(8);
  std::cout << "space time discretization error: "
            << std::setw(8) << gv.size(0) << " elements "
            << std::scientific << l2interpolationerror(u,gfs,x,8) << std::endl;
  {
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U<GV,Real> >(u,"exact solution"));
    vtkwriter.write("instationarytest_exact",Dune::VTK::appendedraw);
  }
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

    if (argc!=3)
      {
        if(helper.rank()==0)
          std::cout << "usage: ./instationarytest <t_level> <x_level>" << std::endl;
        return 1;
      }

    int t_level;
    sscanf(argv[1],"%d",&t_level);

    int x_level;
    sscanf(argv[2],"%d",&x_level);

    // sequential version
    {
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(16));
      std::bitset<dim> periodic(false);
      int overlap=0;
      Dune::YaspGrid<dim> grid(L,N,periodic,overlap);
      grid.globalRefine(x_level);
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      sequential(gv,t_level);
    }
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
