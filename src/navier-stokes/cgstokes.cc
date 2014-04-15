// -*- tab-width: 4; indent-tabs-mode: nil -*-

/** \file

    \brief Example applications of the local operator TaylorHoodNavierStokesJacobian (stationary case).

    This file provides examples applications of Navier-Stokes flow in
    2- and 3-dimensional tubes with and without obstacles. An L-shape
    domain example is provided as well.

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/float_cmp.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/superlu.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/common/parametertreeparser.hh>

#include "../utility/gridexamples.hh"
#include "cgstokes_initial.hh"

//===============================================================
// The driver for all examples
//===============================================================

template<typename GV, typename V_FEM, typename P_FEM, typename IF, typename PRM, int q>
void navierstokes
(
 const GV& gv,
 std::string filename,
 const PRM & parameters,
 V_FEM & vFem, P_FEM & pFem,
 IF & initial_solution )
{
  static const unsigned int dim = GV::dimension;

  typedef double RF;

  Dune::Timer timer;
  std::cout << "=== Initialize:" << timer.elapsed() << std::endl;
  timer.reset();

  ///////////////////////////////////////////////////////
  // Construct grid function spaces
  typedef Dune::PDELab::ConformingDirichletConstraints ConstraintsAssembler;
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::no_blocking,1>
    VectorBackend;
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,V_FEM,dim,VectorBackend,VectorBackend,ConstraintsAssembler> PGFS_V_GFS;
  PGFS_V_GFS powerVGfs(gv,vFem);
  powerVGfs.name("velocity");

  typedef Dune::PDELab::GridFunctionSpace
    <GV, P_FEM, ConstraintsAssembler, VectorBackend> P_GFS;
  P_GFS pGfs(gv,pFem);
  pGfs.name("pressure");

  typedef Dune::PDELab::CompositeGridFunctionSpace
    <VectorBackend,Dune::PDELab::LexicographicOrderingTag, PGFS_V_GFS, P_GFS> GFS;
  GFS gfs(powerVGfs, pGfs);
  ///////////////////////////////////////////////////////

  // Make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();

  // create Taylor-Hood constraints from boundary-type
  typedef Dune::PDELab::StokesVelocityDirichletConstraints<PRM>
    ScalarVelocityConstraints;
  typedef Dune::PDELab::PowerConstraintsParameters<ScalarVelocityConstraints,dim>
    VelocityConstraints;
  typedef Dune::PDELab::StokesPressureDirichletConstraints<PRM>
    PressureConstraints;
  typedef Dune::PDELab::CompositeConstraintsParameters<VelocityConstraints,PressureConstraints>
    Constraints;

  ScalarVelocityConstraints scalarvelocity_constraints(parameters);
  VelocityConstraints velocity_constraints(scalarvelocity_constraints);
  PressureConstraints pressure_constraints(parameters);
  Constraints constraints(velocity_constraints,pressure_constraints);

  Dune::PDELab::constraints(constraints,gfs,cg);

  // Make grid function operator
  typedef Dune::PDELab::TaylorHoodNavierStokes<PRM> LOP;
  LOP lop(parameters,q);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator
    <GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO;
  GO go(gfs,cg,gfs,cg,lop,mbe);

  // Make coefficent vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  Dune::PDELab::interpolate(initial_solution,gfs,x0);
  std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
  timer.reset();

  // Set non constrained dofs to zero
  Dune::PDELab::set_shifted_dofs(cg,0.0,x0);

  // Linear solver
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
  LinearSolver ls(false);

  // Solve (possibly) nonlinear problem
  std::cout << "=== Begin Newton:" << std::endl;
  timer.reset();
  Dune::PDELab::Newton<GO,LinearSolver,V> newton(go,x0,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(30);
  newton.apply();
  std::cout << "=== Finished Newton:" << timer.elapsed() << std::endl;

  // Check residual
  V r(gfs); r=0.;
  go.residual(x0,r);
  std::cout << "Final Residual: " << r.two_norm() << std::endl;


  { // Output grid function with SubsamplingVTKWriter (old school)

    // Generate functions suitable for VTK output
    typedef typename Dune::PDELab::GridFunctionSubSpace
      <GFS,Dune::TypeTree::TreePath<0> > VelocitySubGFS;
    VelocitySubGFS velocitySubGfs(gfs);
    typedef typename Dune::PDELab::GridFunctionSubSpace
      <GFS,Dune::TypeTree::TreePath<1> > PressureSubGFS;
    PressureSubGFS pressureSubGfs(gfs);
    typedef Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,V> VDGF;
    VDGF vdgf(velocitySubGfs,x0);
    typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,V> PDGF;
    PDGF pdgf(pressureSubGfs,x0);

    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
    vtkwriter.write(filename + std::string("_dgf"),Dune::VTK::appendedraw);
  }

  { // Output grid function with SubsamplingVTKWriter

    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x0);
    vtkwriter.write(filename,Dune::VTK::appendedraw);
  }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{

#if !HAVE_SUPERLU
  std::cerr << "Error: These examples work only if SuperLU is available." << std::endl;
  exit(1);
#endif

    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    std::string example_switch;
    if(argc != 2){
      std::cout << std::endl << "PDELab Stokes Examples" << std::endl
                << "----------------------" << std::endl << std::endl
                << "Call with command line parameter to execute examples:" << std::endl << std::endl
                << "              Setup                      | Parameter" << std::endl
                << "-----------------------------------------------------" << std::endl
                << "Hagen-Pouseuille 2D - YaspGrid   - Q2/Q1 :   HY2" << std::endl
                << "Hagen-Pouseuille 2D - AluSimplex - P2/P1 :   HA2" << std::endl
                << "Hagen-Pouseuille 3D - UG - P2/P1         :   HU3" << std::endl
                << "Turbulence Tube  2D - UG - P2/P1         :   TU2" << std::endl
                << "Turbulence Tube  3D - UG - P2/P1         :   TU3" << std::endl
                << "L-Shape Domain   2D - UG - P2/P1         :   LU2" << std::endl
                << std::endl << std::endl
                << "You might also want to take a look at the configuration file \"cgstokes.ini\"."
                << std::endl << std::endl;
      exit(1);
    }
    else
      example_switch = argv[1];

    // Initialize Navier Stokes parameter class from file
    Dune::ParameterTree configuration;

    const std::string config_filename("cgstokes.ini");
    std::cout << "Reading configuration file \""<< config_filename
              << "\"" << std::endl;
    try{
      Dune::ParameterTreeParser::readINITree(config_filename, configuration);
    }
    catch(...){
      std::cerr << "The configuration file \"cgstokes.ini\" "
        "could not be read. Exiting..." << std::endl;
      exit(1);
    }

  try{

    // Yasp Grid Hagen-Poiseuille test 2D
    if(example_switch.find("HY2") != std::string::npos)
    {
      typedef double RF;
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      const int p=2;
      const int q=2*p;

      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,p> V_FEM;
      V_FEM vFem(gv);
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,p-1> P_FEM;
      P_FEM pFem(gv);

      typedef HagenPoiseuilleVelocity<GV,RF,dim> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_HagenPoiseuille BoundaryFunction;
      BoundaryFunction boundary_function;
      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv);
      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"hagenpoiseuille_yasp_Q2Q1_2d", parameters, vFem, pFem, initial_solution);
    }

#if HAVE_ALUGRID
    // ALU Grid Hagen-Poiseuille test 2D
    if(example_switch.find("HA2") != std::string::npos)
    {
      // make grid
      const int dim = 2;
      typedef ALUUnitCube<dim> UnitCube;
      UnitCube unitcube;
      unitcube.grid().globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef UnitCube::GridType::LeafGridView GV;
      const GV& gv=unitcube.grid().leafGridView();

      // make finite element map
      typedef UnitCube::GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;

      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;
      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef double RF;

      typedef HagenPoiseuilleVelocity<GV,RF,dim> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_HagenPoiseuille BoundaryFunction;
      BoundaryFunction boundary_function;
      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv);
      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"hagenpoiseuille_alu_P2P1_2d", parameters, vFem, pFem, initial_solution);
    }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 2D
    if(example_switch.find("TU2") != std::string::npos)
    {
      const int dim = 2;
      typedef Dune::UGGrid<dim> GridType;
      GridType grid;

      typedef double RF;

      std::string grid_file = "grids/turbtube2d.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
      factory.createGrid();

      grid.globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;

      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;
      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef TU_Velocity<GV,RF,2> InitialVelocity;

      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      InitialVelocity init_velocity(gv);
      init_velocity.setTime(1.0); // Take v(t=1.0) as initial velocity!

      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_TU BoundaryFunction;
      typedef PressureDropFlux<GV,RF> NeumannFlux;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 5.0;
      const RF tube_origin = 0.0;

      const RF boundary_pressure = configuration.get<double>("boundaries.pressure");
      BoundaryFunction boundary_function;
      NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);

      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"turbtube_ug_P2P1_2d", parameters, vFem, pFem, initial_solution);
    }
#endif

#if HAVE_UG
    // UG Grid L-shape domain test 2D
    if(example_switch.find("LU2") != std::string::npos)
    {
      const int dim = 2;
      typedef Dune::UGGrid<dim> GridType;
      GridType grid;

      typedef double RF;

      std::string grid_file = "grids/lshape.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
      factory.createGrid();

      grid.globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;

      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;
      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef LU_Velocity<GV,RF,2> InitialVelocity;

      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      InitialVelocity init_velocity(gv);
      init_velocity.setTime(1.0); // Take v(t=1.0) as initial velocity!

      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_LU BoundaryFunction;
      typedef PressureDropFlux<GV,RF> NeumannFlux;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 6.0;
      const RF tube_origin = -1.0;

      const RF boundary_pressure = configuration.get<double>("boundaries.pressure");
      BoundaryFunction boundary_function;
      NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);

      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"lshape_ug_P2P1_2d", parameters, vFem, pFem, initial_solution);
    }
#endif


#if HAVE_UG
    // UG Grid Hagen-Poiseuille test 3D
    if(example_switch.find("HU3") != std::string::npos)
    {
      const int dim = 3;
      typedef Dune::UGGrid<3> GridType;
      GridType grid;

      typedef double RF;

      std::string grid_file = "grids/pipe.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
      factory.createGrid();

      grid.globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;

      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;
      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef HagenPoiseuilleVelocity<GV,RF,3> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_HagenPoiseuille BoundaryFunction;
      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;

      BoundaryFunction boundary_function;
      NeumannFlux neumann_flux(gv);

      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"hagenpoiseuille_ug_P2P1_3d", parameters, vFem, pFem, initial_solution);
    }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 3D
    if(example_switch.find("TU3") != std::string::npos)
    {
      const int dim = 3;
      typedef Dune::UGGrid<dim> GridType;
      GridType grid;

      std::string grid_file = "grids/turbtube.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
      factory.createGrid();
      grid.globalRefine(configuration.get<int>("domain.level"));

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;

      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;
      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef double RF;

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef TU_Velocity<GV,RF,3> InitialVelocity;

      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
        InitialSolution;

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_TU BoundaryFunction;
      typedef PressureDropFlux<GV,RF> NeumannFlux;

      // Domain parameters:
      const int tube_direction = 2; // Tube in z-axes direction
      const RF tube_length = 5.0;
      const RF tube_origin = 0.0;

      const RF boundary_pressure = configuration.get<double>("boundaries.pressure");
      BoundaryFunction boundary_function;
      NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);

      typedef ZeroVectorFunction<GV,RF,dim> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (configuration.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
        (gv,"turbtube_ug_P2P1_3d", parameters, vFem, pFem, initial_solution);

    }
#endif


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
