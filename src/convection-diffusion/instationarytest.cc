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

#include<dune/common/mpihelper.hh>
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
#include<dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/diffusion.hh>
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
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  int
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 1;
    typename Traits::RangeType global = is.geometry().global(x);
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
  //typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,GV::Grid::dimension> FEM;
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<Coord,Real> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
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

  // <<<5>>> Make grid operator space for time-dependent problem
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP; 
  LOP lop(param,4);
  typedef Dune::PDELab::L2 MLOP; 
  MLOP mlop(4);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  //Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::Alexander3Parameter<Real> method;
  typedef typename GFS::template VectorContainer<Real>::Type V;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,MLOP,C,C,MBE> IGOS;
  IGOS igos(method,gfs,cg,gfs,cg,lop,mlop);

  // <<<6>>> Make a linear solver 
// #if HAVE_SUPERLU
//   typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
//   LS ls(false);
// #else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);
  //#endif

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<IGOS,LS,V> PDESOLVER;
  PDESOLVER tnewton(igos,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(0);
  tnewton.setReduction(0.9);
  tnewton.setMinLinearReduction(1e-9);

  // <<<8>>> time-stepper
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,V,V> osm(method,igos,tnewton);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("instationarytest_Q1");
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
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
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
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
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U<GV,Real> >(u,"exact solution"));
    vtkwriter.write("instationarytest_exact",Dune::VTKOptions::binaryappended);
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
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(16);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(x_level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
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
