// -*- tab-width: 4; indent-tabs-mode: nil -*-

/** \file

    \brief Example applications of the local operator TaylorHoodNavierStokesJacobian (instationary case).

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
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/common/parametertreeparser.hh>

#include<dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/stationary/linearproblem.hh>


#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include "../utility/gridexamples.hh"
#include "cgstokes_initial.hh"

//===============================================================
// The driver for all examples
//===============================================================

template<int q, typename GV, typename V_FEM, typename P_FEM, typename IF, typename PRM>
void navierstokes
(
 const GV& gv,
 std::string filename,
 const PRM & parameters,
 const Dune::ParameterTree parser,
 V_FEM & vFem, P_FEM & pFem,
 IF & initial_solution)
{
  static const unsigned int dim = GV::dimension;

  typedef double RF;
  typedef double Real;

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
  LOP lop(parameters);
  typedef Dune::PDELab::NavierStokesMass<PRM> MLOP;
  MLOP mlop(parameters,q);

  Dune::PDELab::FractionalStepParameter<RF> method;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(50); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO0;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);

  typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,Real,Real,Real,C,C> GO1;
  GO1 go1(gfs,cg,gfs,cg,mlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  igo.divideMassTermByDeltaT();
  typedef typename IGO::Traits::Domain V;

  //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // Linear solver
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
  LinearSolver ls(false);

  // Solve (possibly) nonlinear problem
  std::cout << "=== Begin Newton:" << std::endl;
  timer.reset();
  typedef Dune::PDELab::Newton<IGO,LinearSolver,V> PDESolver;
  PDESolver newton(igo,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setMaxIterations(50);
  newton.setLineSearchMaxIterations(30);

  Dune::PDELab::OneStepMethod<Real,IGO,PDESolver,V,V> osm(method,igo,newton);
  osm.setVerbosityLevel(2);

  // Make coefficent vector and initialize it from a function
  V xold(gfs);
  xold = 0.0;
  Dune::PDELab::interpolate(initial_solution,gfs,xold);
  std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
  timer.reset();

  Dune::PDELab::FilenameHelper fn(filename.c_str());
  {
    // Generate functions suitable for VTK output
    typedef typename Dune::PDELab::GridFunctionSubSpace
      <GFS,Dune::TypeTree::TreePath<0> > VelocitySubGFS;
    VelocitySubGFS velocitySubGfs(gfs);
    typedef typename Dune::PDELab::GridFunctionSubSpace
      <GFS,Dune::TypeTree::TreePath<1> > PressureSubGFS;
    PressureSubGFS pressureSubGfs(gfs);
    typedef Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,V> VDGF;
    VDGF vdgf(velocitySubGfs,xold);
    typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,V> PDGF;
    PDGF pdgf(pressureSubGfs,xold);

    // Output grid function with SubsamplingVTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
    fn.increment();
  }

  timer.reset();

  Real time = 0;
  Real final_time = parser.get("temporal.time",double(1.0));
  Real dt = parser.get("temporal.tau",double(0.1));
  Real dt_min = 1e-6;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(initial_solution,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);
  while (time < final_time - dt_min*0.5)
    {
      // do time step
      osm.apply(time,dt,xold,initial_solution,x);

      {
        // Generate functions suitable for VTK output
        typedef typename Dune::PDELab::GridFunctionSubSpace
          <GFS,Dune::TypeTree::TreePath<0> > VelocitySubGFS;
        VelocitySubGFS velocitySubGfs(gfs);
        typedef typename Dune::PDELab::GridFunctionSubSpace
          <GFS,Dune::TypeTree::TreePath<1> > PressureSubGFS;
        PressureSubGFS pressureSubGfs(gfs);
        typedef Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,V> VDGF;
        VDGF vdgf(velocitySubGfs,x);
        typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,V> PDGF;
        PDGF pdgf(pressureSubGfs,x);

        // Output grid function with SubsamplingVTKWriter
        Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
        vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
        fn.increment();
      }

      xold = x;
      time += dt;
    }

  std::cout << "=== Total time:" << timer.elapsed() << std::endl;

  // Compute norm of final solution
  typedef typename Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<0> > VelocitySubGFS;
  VelocitySubGFS velocitySubGfs(gfs);
  typedef typename Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<1> > PressureSubGFS;
  PressureSubGFS pressureSubGfs(gfs);
  typedef Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,V> VDGF;
  VDGF vdgf(velocitySubGfs,x);
  typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,V> PDGF;
  PDGF pdgf(pressureSubGfs,x);
  typename PDGF::Traits::RangeType l1norm(0);
  Dune::PDELab::integrateGridFunction(pdgf,l1norm,q);
  std::cout << std::setw(12) << std::setprecision(7) << std::scientific
            << "L1 norm: " <<  l1norm << std::endl;


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
                << "Turbulence Tube  2D - UG - P2/P1         :   TU2" << std::endl
                << "L-Shape Domain   2D - UG - P2/P1         :   LU2" << std::endl
                << std::endl << std::endl
                << "You might also want to take a look at the configuration file \"cgstokes_instat.ini\"."
                << std::endl << std::endl;
      exit(1);
    }
    else
      example_switch = argv[1];

    // Initialize Navier Stokes parameter class from file
    Dune::ParameterTree config_parser;
    const std::string config_filename("cgstokes_instat.ini");
    std::cout << "Reading configuration file \""<< config_filename
              << "\"" << std::endl;
    try{
      Dune::ParameterTreeParser::readINITree(config_filename,config_parser);
    }
    catch(...){
      std::cerr << "The configuration file \"cgstokes_instat.ini\" "
        "could not be read. Exiting..." << std::endl;
      exit(1);
    }

  try{

#if HAVE_UG
    // UG Grid turbulence tube test 2D
    if(example_switch.find("TU2") != std::string::npos)
    {
      typedef Dune::UGGrid<2> GridType;
      GridType grid;

      typedef double RF;

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/turbtube2d.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,boundary_index_map,element_index_map,true,false);
      factory.createGrid();
      grid.globalRefine(config_parser.get<int>("domain.level"));

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

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(gv);

      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_TU BoundaryFunction;
      typedef PressureDropFlux<GV,RF> NeumannFlux;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 5.0;
      const RF tube_origin = 0.0;

      BoundaryFunction boundary_function;
      const RF boundary_pressure = config_parser.get<double>("boundaries.pressure");
      NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);
      typedef ZeroVectorFunction<GV,RF,2> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (config_parser.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<q>
        (gv,"turbtube_ug_P2P1_2d", parameters, config_parser, vFem, pFem, initial_solution);
    }
#endif

#if HAVE_UG
    // UG Grid L-shape domain test 2D
    if(example_switch.find("LU2") != std::string::npos)
    {
      typedef Dune::UGGrid<2> GridType;
      GridType grid;

      typedef double RF;

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/lshape.msh";
      Dune::GridFactory<GridType> factory(&grid);
      Dune::GmshReader<GridType>::read(factory,grid_file,boundary_index_map,element_index_map,true,false);
      factory.createGrid();
      grid.globalRefine(config_parser.get<int>("domain.level"));

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

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef BCTypeParam_TU BoundaryFunction;
      typedef PressureDropFlux<GV,RF> NeumannFlux;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 6.0;
      const RF tube_origin = -1.0;

      BoundaryFunction boundary_function;
      const RF boundary_pressure = config_parser.get<double>("boundaries.pressure");
      NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);
      typedef ZeroVectorFunction<GV,RF,2> SourceFunction;
      SourceFunction source_function(gv);

      typedef Dune::PDELab::NavierStokesDefaultParameters
        <GV,RF,SourceFunction,BoundaryFunction,InitialSolution,NeumannFlux>
        LOPParameters;
      LOPParameters parameters
        (config_parser.sub("physics"),source_function,boundary_function,
         initial_solution,neumann_flux);

      // solve problem
      navierstokes<q>
        (gv,"lshape_ug_P2P1_2d", parameters, config_parser, vFem, pFem, initial_solution);
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
