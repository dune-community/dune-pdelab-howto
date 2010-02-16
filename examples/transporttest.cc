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
#include<dune/pdelab/finiteelementmap/p0constraints.hh>
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
#include<dune/pdelab/localoperator/transportccfv.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include"gridexamples.hh"

//==============================================================================
// Parameter class for the convection diffusion problem
//==============================================================================

//! base class for parameter class
template<typename GV, typename RF>
class TransportProblem : 
  public Dune::PDELab::TransportSpatialParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>, 
                                                          TransportProblem<GV,RF> >
{
public:
  typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType flux;
    flux[0] = 1.0;
    flux[1] = 0.5;
    return flux;
  }

  //! tensor permeability
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1E-10;
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
    if (global[0]<1E-6 || global[1]<1E-6)
      return 1;
    else return 2;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType global = e.geometry().global(x);
    if (global[0]>1E-12 && global[0]<1.0-1E-12 &&           // this is needed for parallel
        global[1]>1E-12 && global[1]<1.0-1E-12) return 0.0; // overlapping version !
    if (global[0]<1E-6 && global[1]>0.25 && global[1]<0.75)
      return 1.0;
    else
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

  //! capacity function
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 1.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {}
};


//===============================================================
// Some variants to solve the nonlinear diffusion problem
//===============================================================

// example solving stationary transport problem
template<class GV>
void stationary (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(Dune::GeometryType::cube); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_Transport<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_Transport<Param> G;
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
  typedef Dune::PDELab::CCFVSpatialTransportOperator<Param> LOP; 
  LOP lop(param);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
#ifdef HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls(false);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(0.9);
  newton.setMinLinearReduction(1e-9);
  newton.apply();

  // <<<8>>> graphical output
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,x);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write("transporttest_stationary",Dune::VTKOptions::binaryappended);
  }
}


// example using implicit time-stepping
template<class GV>
void implicit_scheme (const GV& gv, double Tend, double timestep)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(Dune::GeometryType::cube); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_Transport<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_Transport<Param> G;
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
  typedef Dune::PDELab::CCFVSpatialTransportOperator<Param> LOP; 
  LOP lop(param);
  typedef Dune::PDELab::CCFVTemporalOperator<Param> SLOP; 
  SLOP slop(param);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);
  Dune::PDELab::FractionalStepParameter<Real> method;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,SLOP,C,C,MBE> IGOS;
  IGOS igos(method,gfs,cg,gfs,cg,lop,slop);

  // <<<6>>> Make a linear solver
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,C> LS;
  LS ls(gfs,cg,5000,1,1); 

  // <<<7>>> make Newton for time-dependent problem
  typedef Dune::PDELab::Newton<IGOS,LS,V> PDESOLVER;
  PDESOLVER tnewton(igos,ls);
  tnewton.setReassembleThreshold(0.0);
  tnewton.setVerbosityLevel(0);
  tnewton.setReduction(0.9);
  tnewton.setMinLinearReduction(1e-7);

  // <<<8>>> time-stepper
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,V,V> osm(method,igos,tnewton);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("transporttest_implicit_fs");
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  x = 0;
  while (time < Tend)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,x);
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
      vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
      fn.increment();

      xold = x;
      time += dt;
    }
}


// example using explicit time-stepping
template<class GV>
void explicit_scheme (const GV& gv, double Tend, double timestep)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(Dune::GeometryType::cube); // works only for cubes
  typedef Dune::PDELab::P0ParallelConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef TransportProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_Transport<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_Transport<Param> G;
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
  typedef Dune::PDELab::CCFVSpatialTransportOperator<Param> LOP; 
  LOP lop(param);
  typedef Dune::PDELab::CCFVTemporalOperator<Param> SLOP; 
  SLOP slop(param);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);
  Dune::PDELab::ExplicitEulerParameter<Real> method;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,SLOP,C,C,MBE> IGOS;
  IGOS igos(method,gfs,cg,gfs,cg,lop,slop);

  // <<<6>>> Make a linear solver backend
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<8>>> time-stepper
  typedef Dune::PDELab::CFLTimeController<Real,IGOS> TC;
  TC tc(0.999,igos);
  Dune::PDELab::ExplicitOneStepMethod<Real,IGOS,LS,V,V,TC> osm(method,igos,ls,tc);
  osm.setVerbosityLevel(2);

  // <<<9>>> initial value and initial value for first time step with b.c. set
  V xold(gfs,0.0);
  xold = 0.0;

  // <<<10>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("transporttest_explicit_euler");
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  x = 0;
  while (time < Tend)
    {
      // do time step
      dt = osm.apply(time,dt,xold,x);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,x);
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
      vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
      fn.increment();

      xold = x;
      time += dt;
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

	if (argc!=5)
	  {
		if(helper.rank()==0)
		  std::cout << "usage: ./transporttest <end time> <time step> <elements on a side> <overlap>" << std::endl;
		return 1;
	  }

	double Tend;
	sscanf(argv[1],"%lg",&Tend);

	double timestep;
	sscanf(argv[2],"%lg",&timestep);

	int n;
	sscanf(argv[3],"%d",&n);

	int o;
	sscanf(argv[4],"%d",&o);

#if HAVE_MPI
    // parallel overlapping version
    if (1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(n);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=o;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      //stationary(gv);
      //implicit_scheme(gv,Tend,timestep);
      explicit_scheme(gv,Tend,timestep);
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
