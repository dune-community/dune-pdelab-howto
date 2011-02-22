// -*- tab-width: 4; indent-tabs-mode: nil -*-

/** \file 

    \brief Example applications of the local operator
    TaylorHoodNavierStokesJacobian. 

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
#include<dune/common/mpihelper.hh>
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

#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/finiteelementmap/hangingnodeconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/common/configparser.hh>

#include "../utility/gridexamples.hh"
#include "cg_stokes_initial.hh"

//! A simple parameter class which complies with the interface
//! required by the local operator. It is initialized with a
//! Dune::ConfigParser object.
class NavierStokesParameters 
  : public Dune::PDELab::TaylorHoodNavierStokesParameters<double>
{
public:  
  double rho_;
  double mu_;
  double pressure;
  int domain_level;

  double rho() const{
    return rho_;
  }

  double mu() const{
    return mu_;
  }

  template <class ConfigParser>
  NavierStokesParameters(ConfigParser & parser)
    : rho_(1.0), mu_(1.0), pressure(1.0), domain_level(1)
  {
    bool valid = true;
    valid = valid & parser.hasKey("physics.mu");
    valid = valid & parser.hasKey("physics.rho");
    valid = valid & parser.hasKey("boundaries.pressure");
    valid = valid & parser.hasKey("domain.level");

    if(!valid){
      std::cerr << "Error: The configuration file "
                << "was not valid. Default values will be used."
                << std::endl;
    }
    else{
      rho_ = parser.get("physics.rho",double(1.0));
      mu_ = parser.get("physics.mu",double(1.0));
      pressure = parser.get("boundaries.pressure",double(1.0));
      domain_level = parser.get("domain.level",int(1));
    }
  }
};

//===============================================================
// The driver for all examples
//===============================================================

template<typename GV, typename V_FEM, typename P_FEM, typename IF, typename BF, typename NF, int q> 
void navierstokes 
(
 const GV& gv, 
 std::string filename, 
 const NavierStokesParameters & parameters,
 V_FEM & vFem, P_FEM & pFem, 
 IF & initial_solution,
 BF & boundary_function,
 NF & neumann_flux)
{
  typedef typename GV::Grid::ctype DF;
  static const unsigned int dim = GV::dimension;

  typedef double RF;
  typedef BF BoundaryFunction;
  typedef IF InitializationFunction;
  typedef NF NeumannFlux;

  Dune::Timer timer;
  std::cout << "=== Initialize:" << timer.elapsed() << std::endl;
  timer.reset();

  ///////////////////////////////////////////////////////
  // Construct grid function spaces
  typedef Dune::PDELab::ISTLVectorBackend<1> VectorBackend;
  typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
  typedef Dune::PDELab::SimpleGridFunctionStaticSize GFSSize;
  typedef Dune::PDELab::GridFunctionSpace
    <GV, V_FEM, Constraints, VectorBackend, GFSSize> V_GFS;
  V_GFS vGfs(gv,vFem);

  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;
  typedef Dune::PDELab::PowerGridFunctionSpace<V_GFS,dim,GFMapper> PGFS_V_GFS;
  PGFS_V_GFS powerVGfs(vGfs);

  typedef Dune::PDELab::GridFunctionSpace
    <GV, P_FEM, Constraints, VectorBackend, GFSSize> P_GFS;
  P_GFS pGfs(gv,pFem);

  typedef Dune::PDELab::CompositeGridFunctionSpace<GFMapper,PGFS_V_GFS, P_GFS> GFS;
  GFS gfs(powerVGfs, pGfs);
  ///////////////////////////////////////////////////////

  // Make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(boundary_function,gfs,cg);

  // Make coefficent vector and initialize it from a function
  typedef typename GFS::template VectorContainer<RF>::Type V;
  V x0(gfs);
  x0 = 0.0;
  Dune::PDELab::interpolate(initial_solution,gfs,x0);
  std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
  timer.reset();

  // Set non constrained dofs to zero
  Dune::PDELab::set_shifted_dofs(cg,0.0,x0);

  // Make grid function operator
  typedef Dune::PDELab::TaylorHoodNavierStokesJacobian
    <BoundaryFunction,NeumannFlux,NavierStokesParameters,true,q> 
    LOP; 
  LOP lop(boundary_function,neumann_flux,parameters);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // Linear solver
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
  LinearSolver ls(false);

  // Solve (possibly) nonlinear problem
  std::cout << "=== Begin Newton:" << std::endl;
  timer.reset();
  Dune::PDELab::Newton<GOS,LinearSolver,V> newton(gos,x0,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(30);
  newton.apply();
  std::cout << "=== Finished Newton:" << timer.elapsed() << std::endl;

  // Check residual
  V r(gfs); r=0.;
  gos.residual(x0,r);
  std::cout << "Final Residual: " << r.two_norm() << std::endl;

  // Generate functions suitable for VTK output
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> VelocitySubGFS;
  VelocitySubGFS velocitySubGfs(gfs);
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,1> PressureSubGFS;
  PressureSubGFS pressureSubGfs(gfs);
  typedef Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,V> VDGF;
  VDGF vdgf(velocitySubGfs,x0);
  typedef Dune::PDELab::DiscreteGridFunction<PressureSubGFS,V> PDGF;
  PDGF pdgf(pressureSubGfs,x0);

  // Output grid function with SubsamplingVTKWriter
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
  vtkwriter.write(filename,Dune::VTKOptions::binaryappended);
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
                << "You might also want to take a look at the configuration file \"cg_stokes.ini\"." 
                << std::endl << std::endl;
      exit(1);
    }
    else
      example_switch = argv[1];

    // Initialize Navier Stokes parameter class from file
    Dune::ConfigParser config_parser;
    const std::string config_filename("cg_stokes.ini");
    std::cout << "Reading configuration file \""<< config_filename 
              << "\"" << std::endl;
    try{
      config_parser.parseFile(config_filename,false);
    }
    catch(...){
      std::cerr << "The configuration file \"cg_stokes.ini\" "
        "could not be read. Exiting..." << std::endl;
      exit(1);
    }
    NavierStokesParameters parameters(config_parser);

  try{

    // Yasp Grid Hagen-Poiseuille test 2D
    if(example_switch.find("HY2") != std::string::npos)
    {
      typedef double RF;
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(parameters.domain_level);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView(); 

      const int p=2;
      const int q=2*p;

      typedef Dune::PDELab::Q22DLocalFiniteElementMap<double,double> V_FEM;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<double,double,2> P_FEM;

      V_FEM vFem; P_FEM pFem;

      typedef HagenPoiseuilleVelocity<GV,RF,2> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef HagenPoiseuilleVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,2> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv);

      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"hagenpoiseuille_yasp_Q2Q1_2d", parameters, vFem, pFem, initial_solution, 
         boundary_function, neumann_flux);
    }

#if HAVE_ALUGRID
    // ALU Grid Hagen-Poiseuille test 2D
    if(example_switch.find("HA2") != std::string::npos)
    {
      // make grid
      typedef ALUUnitCube<2> UnitCube;
      UnitCube unitcube;
      unitcube.grid().globalRefine(parameters.domain_level);

      // get view
      typedef UnitCube::GridType::LeafGridView GV;
      const GV& gv=unitcube.grid().leafView(); 
 
      // make finite element map
      typedef UnitCube::GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;

      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef double RF;

      typedef HagenPoiseuilleVelocity<GV,RF,2> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef HagenPoiseuilleVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,2> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv);
  
      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"hagenpoiseuille_alu_P2P1_2d", parameters, vFem, pFem, initial_solution, 
         boundary_function, neumann_flux);
    }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 2D
    if(example_switch.find("TU2") != std::string::npos)
    {
      typedef Dune::UGGrid<2> GridType;
      GridType grid(400);

      typedef double RF;

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/turbtube2d.msh";
      Dune::GmshReader<GridType> gmsh_reader;
      gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);
      grid.globalRefine(parameters.domain_level);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;

      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef Dune::PDELab::PowerGridFunction<ZeroFunction,2> 
        InitialVelocity;
      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(zero_function);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef PressureDropVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,2> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 5.0;
      const RF tube_origin = 0.0;


      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv, tube_length, tube_origin, tube_direction);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef PressureDropFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv, parameters.pressure, tube_length, tube_origin, tube_direction);
  
      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"turbtube_ug_P2P1_2d", parameters, vFem, pFem, 
         initial_solution, boundary_function,neumann_flux);
    }
#endif

#if HAVE_UG
    // UG Grid L-shape domain test 2D
    if(example_switch.find("LU2") != std::string::npos)
    {
      typedef Dune::UGGrid<2> GridType;
      GridType grid(400);

      typedef double RF;

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/lshape.msh";
      Dune::GmshReader<GridType> gmsh_reader;
      gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);
      grid.globalRefine(parameters.domain_level);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;

      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef Dune::PDELab::PowerGridFunction<ZeroFunction,2> 
        InitialVelocity;
      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(zero_function);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef PressureDropVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,2> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      // Domain parameters:
      const int tube_direction = 0; // Tube in x-axes direction
      const RF tube_length = 6.0;
      const RF tube_origin = -1.0;


      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv, tube_length, tube_origin, tube_direction);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef PressureDropFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv, parameters.pressure, tube_length, tube_origin, tube_direction);
  
      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"lshape_ug_P2P1_2d", parameters, vFem, pFem, 
         initial_solution, boundary_function,neumann_flux);
    }
#endif


#if HAVE_UG
    // UG Grid Hagen-Poiseuille test 3D
    if(example_switch.find("HU3") != std::string::npos)
    {
      typedef Dune::UGGrid<3> GridType;
      GridType grid(400);

      typedef double RF;

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/pipe.msh";
      Dune::GmshReader<GridType> gmsh_reader;
      gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);
      grid.globalRefine(parameters.domain_level);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;

      V_FEM vFem(gv); P_FEM pFem(gv);

      typedef HagenPoiseuilleVelocity<GV,RF,3> InitialVelocity;
      typedef ZeroScalarFunction<GV,RF> InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      InitialVelocity init_velocity(gv);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);

      typedef HagenPoiseuilleVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,3> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv);
  
      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"hagenpoiseuille_ug_P2P1_3d", parameters, vFem, pFem, 
         initial_solution, boundary_function,neumann_flux);
    }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 3D
    if(example_switch.find("TU3") != std::string::npos)
    {
      typedef Dune::UGGrid<3> GridType;
      GridType grid(400);

      std::vector<int> boundary_index_map;
      std::vector<int> element_index_map;

      std::string grid_file = "grids/turbtube.msh";
      Dune::GmshReader<GridType> gmsh_reader;
      gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);
      grid.globalRefine(parameters.domain_level);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=2;
      const int q=2*k;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k> V_FEM;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k-1> P_FEM;

      V_FEM vFem(gv); P_FEM pFem(gv);
  
      typedef double RF;

      typedef ZeroScalarFunction<GV,RF> ZeroFunction;
      typedef Dune::PDELab::PowerGridFunction<ZeroFunction,3> 
        InitialVelocity;
      typedef ZeroFunction InitialPressure;
      typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure> 
        InitialSolution;

      ZeroFunction zero_function(gv);
      InitialVelocity init_velocity(zero_function);
      InitialPressure init_pressure(gv);
      InitialSolution initial_solution(init_velocity,init_pressure);


      typedef PressureDropVelocityBoundaryFunction<GV> ScalarVelocityBoundaryFunction;
      typedef Dune::PDELab::PowerGridFunction<ScalarVelocityBoundaryFunction,3> 
        VelocityBoundaryFunction;
      typedef ScalarNeumannBoundaryFunction<GV> PressureBoundaryFunction;
      typedef Dune::PDELab::CompositeGridFunction<VelocityBoundaryFunction,PressureBoundaryFunction> 
        BoundaryFunction;

      // Domain parameters:
      const int tube_direction = 2; // Tube in z-axes direction
      const RF tube_length = 5.0;
      const RF tube_origin = 0.0;

      ScalarVelocityBoundaryFunction bf_scalar_velocity(gv, tube_length, tube_origin, tube_direction);
      VelocityBoundaryFunction bf_velocity(bf_scalar_velocity);
      PressureBoundaryFunction bf_pressure(gv);
      BoundaryFunction boundary_function(bf_velocity,bf_pressure);

      typedef PressureDropFlux<GV,RF> NeumannFlux;
      NeumannFlux neumann_flux(gv, parameters.pressure, tube_length, tube_origin, tube_direction);
  
      // solve problem
      navierstokes<GV,V_FEM,P_FEM,InitialSolution,BoundaryFunction,NeumannFlux,q>
        (gv,"turbtube_ug_P2P1_3d", parameters, vFem, pFem, initial_solution, boundary_function, neumann_flux);

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
