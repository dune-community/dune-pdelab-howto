// -*- tab-width: 4; indent-tabs-mode: nil -*-

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

#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>

#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/seqistlsolverbackend.hh>
//#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
//#include <dune/pdelab/gridfunctionspace/entityblockedlocalordering.hh>
//#include <dune/pdelab/gridfunctionspace/leaflocalordering.hh>

#include <dune/pdelab/newton/newton.hh>

#include "../utility/gridexamples.hh"
#include "../navier-stokes/cgstokes_initial.hh"
#include <dune/pdelab/common/benchmarkhelper.hh>

//===============================================================
// The driver for all examples
//===============================================================

template<typename GV, typename V_FEM, typename P_FEM, typename IF, typename PRM, int q>
void navierstokes(const GV& gv,
                  std::string filename,
                  const PRM & parameters,
                  V_FEM & vFem, P_FEM & pFem,
                  IF & initial_solution,
                  bool solve,
                  bool io,
                  std::size_t runs)
{

  Dune::PDELab::BenchmarkHelper<> bh(filename,runs);

  for (std::size_t run = 0; run < runs; ++run)
    {
      bh.start_run(std::cout);

      typedef typename GV::Grid::ctype DF;
      static const unsigned int dim = GV::dimension;

      typedef double RF;
      typedef IF InitializationFunction;

      bh.start("global setup",std::cout);
      bh.start("GFS setup",std::cout);

      ///////////////////////////////////////////////////////
      // Construct grid function spaces
      typedef Dune::PDELab::ISTLVectorBackend<> VectorBackend;
      typedef Dune::PDELab::ISTLVectorBackend<> VelocityVectorBackend;

#ifdef CGSTOKES_MACRO_BLOCKS
      typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::dynamic_blocking> TaylorHoodVectorBackend;
#else
      typedef Dune::PDELab::ISTLVectorBackend<> TaylorHoodVectorBackend;
#endif

#if CGSTOKES_BLOCK_LEVEL > 0
      typedef Dune::PDELab::EntityBlockedOrderingTag VelocityOrderingTag;
#else
      typedef Dune::PDELab::LexicographicOrderingTag VelocityOrderingTag;
#endif

#if CGSTOKES_BLOCK_LEVEL == 2

#ifdef CGSTOKES_MACRO_BLOCKS
#error "You cannot use macro blocks when aggregating the whole function space per entity"
#endif

      typedef Dune::PDELab::EntityBlockedOrderingTag TaylorHoodOrderingTag;
#else
      typedef Dune::PDELab::LexicographicOrderingTag TaylorHoodOrderingTag;
#endif

      typedef Dune::PDELab::ConformingDirichletConstraints ConstraintsAssembler;

      typedef Dune::PDELab::VectorGridFunctionSpace<
        GV,V_FEM,dim,
        VelocityVectorBackend,
        VectorBackend,
        ConstraintsAssembler,
        VelocityOrderingTag> PGFS_V_GFS;
      PGFS_V_GFS powerVGfs(gv,vFem);
      powerVGfs.name("v");

      typedef Dune::PDELab::GridFunctionSpace<
        GV, P_FEM,
        ConstraintsAssembler,
        VectorBackend> P_GFS;
      P_GFS pGfs(gv,pFem);
      pGfs.name("p");

      typedef Dune::PDELab::CompositeGridFunctionSpace<
        TaylorHoodVectorBackend,
        TaylorHoodOrderingTag,
        PGFS_V_GFS,
        P_GFS> GFS;
      GFS gfs(powerVGfs, pGfs);

      bh.end("GFS setup",std::cout);
      bh.start("ordering update",std::cout);

      gfs.ordering();

      bh.end("ordering update",std::cout);
      bh.end("global setup",std::cout);

      ///////////////////////////////////////////////////////

      bh.start("constraints",std::cout);

      // Make constraints map and initialize it from a function
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      C cg;
      cg.clear();

      // create Taylor-Hood constraints from boundary-type
      typedef Dune::PDELab::StokesVelocityDirichletConstraints<PRM>
        ScalarVelocityConstraints;
      typedef Dune::PDELab::PowerConstraintsParameters<ScalarVelocityConstraints,dim>
        VelocityConstraints;
      typedef Dune::PDELab::NoDirichletConstraintsParameters
        PressureConstraints;
      typedef Dune::PDELab::CompositeConstraintsParameters<VelocityConstraints,PressureConstraints>
        Constraints;

      ScalarVelocityConstraints scalarvelocity_constraints(parameters);
      VelocityConstraints velocity_constraints(scalarvelocity_constraints);
      PressureConstraints pressure_constraints;
      Constraints constraints(velocity_constraints,pressure_constraints);

      Dune::PDELab::constraints(constraints,gfs,cg);

      bh.end("constraints",std::cout);
      bh.start("LOP construction",std::cout);

      // Make grid function operator
      typedef Dune::PDELab::TaylorHoodNavierStokes<PRM> LOP;
      LOP lop(parameters,q);

      bh.end("LOP construction",std::cout);
      bh.start("GOP construction",std::cout);

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
      typedef Dune::PDELab::GridOperator
        <GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO;
      GO go(gfs,cg,gfs,cg,lop,mbe);

      bh.end("GOP construction",std::cout);
      bh.start("vector creation",std::cout);

      // Make coefficent vector and initialize it from a function
      typedef typename GO::Traits::Domain V;
      V x0(gfs);
      x0 = 0.0;

      bh.end("vector creation",std::cout);
      bh.start("interpolation",std::cout);

      Dune::PDELab::interpolate(initial_solution,gfs,x0);

      bh.end("interpolation",std::cout);
      bh.start("set shifted DOFs",std::cout);

      // Set non constrained dofs to zero
      Dune::PDELab::set_shifted_dofs(cg,0.0,x0);

      bh.end("set shifted DOFs",std::cout);

      {
        bh.start("matrix creation",std::cout);

        typename GO::Traits::Jacobian A(go);

        bh.end("matrix creation",std::cout);
        bh.start("jacobian",std::cout);

        go.jacobian(x0,A);

        bh.end("jacobian",std::cout);
      }

      bh.start("solve",std::cout);

      if (solve)
        {
          // Linear solver
          typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
          LinearSolver ls(false);

          // Solve (possibly) nonlinear problem
          Dune::PDELab::Newton<GO,LinearSolver,V> newton(go,x0,ls);
          newton.setReassembleThreshold(0.0);
          newton.setVerbosityLevel(2);
          newton.setMaxIterations(25);
          newton.setLineSearchMaxIterations(30);
          newton.apply();
        }

      bh.end("solve",std::cout);

      // Check residual
      V r(gfs); r=0.;

      bh.start("residual",std::cout);

      go.residual(x0,r);

      bh.end("residual",std::cout);

      bh.start("IO",std::cout);

      if (io)
        {
          // Output grid function with SubsamplingVTKWriter
          Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
          Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x0);
          vtkwriter.write(filename,Dune::VTK::appendedraw);
        }

      bh.end("IO",std::cout);

      bh.end_run(std::cout);

    }

  bh.print(std::cout);

  std::stringstream timings_name;
  timings_name << filename << "_timings.txt";

  std::ofstream timings_file(timings_name.str());
  timings_file << filename << " " << runs << std::endl;
  bh.print(timings_file);
}


void merge_sub_trees(const Dune::ParameterTree& source, Dune::ParameterTree& target_father, std::string sub_name)
{
  // Make sure the sub exists within target
  // As we cannot create a sub directly, we just set a dummy key within the sub
  if (!target_father.hasSub(sub_name))
    target_father[sub_name + "._dummy"] = "_dummy";

  Dune::ParameterTree& target = target_father.sub(sub_name);

  typedef Dune::ParameterTree::KeyVector::const_iterator key_iterator_type;

  // copy missing keys
  for (key_iterator_type it = source.getValueKeys().begin(), end = source.getValueKeys().end(); it != end; ++it)
    if (!target.hasKey(*it))
      target[*it] = source[*it];

  // recursively copy missing subs
  for (key_iterator_type it = source.getSubKeys().begin(), end = source.getSubKeys().end(); it != end; ++it)
    if (!target.hasSub(*it))
      merge_sub_trees(source.sub(*it),target,*it);
}


template<typename Grid>
std::string name(std::string problem, const Grid& grid, std::size_t level, bool solve, bool io)
{
  std::stringstream n;
  n << "cgstokes_" << problem
    << "_" << Grid::dimension << "D"
    << "_l" << level
    << (solve ? "_solve" : "_nosolve")
    << (io ? "_io" : "_noio");

#ifdef CGSTOKES_MACROBLOCKS
  n << "_macroblocks";
#endif

#if CGSTOKES_BLOCK_LEVEL == 1
  n << "_velocityblocked";
#elif CGSTOKES_BLOCK_LEVEL == 2
  n << "_allblocked";
#endif

  return n.str();
}


//===============================================================
// Main program with grid setup
//===============================================================


int main(int argc, char** argv)
{

  try {

    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    Dune::ParameterTree params;

    if (argc == 1)
      {
        std::cerr << "no parameter file passed, defaulting to cgstokes.ini..." << std::endl;
        Dune::ParameterTreeParser::readINITree("cgstokes.ini",params);
      }
    else if (argc == 2)
      {
        std::cerr << "reading parameters from " << argv[1] << "..." << std::endl;
        Dune::ParameterTreeParser::readINITree(argv[1],params);
      }
    else
      {
        std::cerr << "Usage: " << argv[0] << " [parameter file]" << std::endl;
        return 64;
      }

    const std::size_t global_runs = params.get("global.runs",5);
    const bool global_solve = params.get("global.solve",false);
    const bool global_io = params.get("global.io",false);
    std::string grid_base_dir = params["global.griddirectory"];

    // Yasp Grid Hagen-Poiseuille test 2D
    if (params.hasSub("HY2") && params.get<bool>("HY2.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("HY2");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");

        typedef double RF;
        // make grid
        const int dim = 2;
        Dune::FieldVector<double,dim> L(1.0);
        Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
        std::bitset<dim> B(false);
        Dune::YaspGrid<dim> grid(L,N,B,0);
        grid.globalRefine(level);

        // get view
        typedef Dune::YaspGrid<dim>::LeafGridView GV;
        const GV& gv=grid.leafGridView();

        typedef Dune::YaspGrid<dim>::ctype DF;

        const int p=2;
        const int q=2*p;

        typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,p> V_FEM;
        V_FEM vFem(gv);
        typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,p-1> P_FEM;
        P_FEM pFem(gv);

        typedef HagenPoiseuilleVelocity<GV,RF,2> InitialVelocity;
        typedef ZeroScalarFunction<GV,RF> InitialPressure;
        typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
          InitialSolution;

        InitialVelocity init_velocity(gv);
        InitialPressure init_pressure(gv);
        InitialSolution initial_solution(init_velocity,init_pressure);

        typedef BCTypeParam_HagenPoiseuille BoundaryFunction;
        typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
        typedef ZeroVectorFunction<GV,RF,2> ExternalForce;

        BoundaryFunction boundary_function;
        NeumannFlux neumann_flux(gv);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("HY2",grid,level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
      }

#if HAVE_ALUGRID
    // ALU Grid Hagen-Poiseuille test 2D
    if (params.hasSub("HA2") && params.get<bool>("HA2.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("HA2");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");

        // make grid
        typedef ALUUnitCube<2> UnitCube;
        UnitCube unitcube;
        unitcube.grid().globalRefine(level);

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

        typedef HagenPoiseuilleVelocity<GV,RF,2> InitialVelocity;
        typedef ZeroScalarFunction<GV,RF> InitialPressure;
        typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
          InitialSolution;

        InitialVelocity init_velocity(gv);
        InitialPressure init_pressure(gv);
        InitialSolution initial_solution(init_velocity,init_pressure);

        typedef BCTypeParam_HagenPoiseuille BoundaryFunction;
        typedef HagenPoiseuilleZeroFlux<GV,RF> NeumannFlux;
        typedef ZeroVectorFunction<GV,RF,2> ExternalForce;

        BoundaryFunction boundary_function;
        NeumannFlux neumann_flux(gv);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("HA2",unitcube.grid(),level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
      }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 2D
    if (params.hasSub("TU2") && params.get<bool>("TU2.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("TU2");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");
        merge_sub_trees(params.sub("boundaries"),lp,"boundaries");

        typedef Dune::UGGrid<2> GridType;
        GridType grid;

        typedef double RF;

        std::string grid_file = grid_base_dir + "/turbtube2d.msh";
        Dune::GridFactory<GridType> factory(&grid);
        Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
        factory.createGrid();

        grid.globalRefine(level);

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
        typedef Dune::PDELab::PowerGridFunction<ZeroFunction,2>
          InitialVelocity;
        typedef ZeroFunction InitialPressure;
        typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
          InitialSolution;

        ZeroFunction zero_function(gv);
        InitialVelocity init_velocity(zero_function);
        InitialPressure init_pressure(gv);
        InitialSolution initial_solution(init_velocity,init_pressure);

        typedef BCTypeParam_PressureDrop<GV> BoundaryFunction;
        typedef PressureDropFlux<GV,RF> NeumannFlux;
        typedef ZeroVectorFunction<GV,RF,2> ExternalForce;

        // Domain parameters:
        const int tube_direction = 0; // Tube in x-axes direction
        const RF tube_length = 5.0;
        const RF tube_origin = 0.0;

        const RF boundary_pressure = lp.get<double>("boundaries.pressure");
        BoundaryFunction boundary_function(tube_length, tube_origin, tube_direction);
        NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("TU2",grid,level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
      }
#endif

#if HAVE_UG
    // UG Grid L-shape domain test 2D
    if (params.hasSub("LU2") && params.get<bool>("LU2.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("LU2");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");
        merge_sub_trees(params.sub("boundaries"),lp,"boundaries");

        typedef Dune::UGGrid<2> GridType;
        GridType grid;

        typedef double RF;

        std::string grid_file = grid_base_dir + "/lshape.msh";
        Dune::GridFactory<GridType> factory(&grid);
        Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
        factory.createGrid();

        grid.globalRefine(level);

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
        typedef Dune::PDELab::PowerGridFunction<ZeroFunction,2>
          InitialVelocity;
        typedef ZeroFunction InitialPressure;
        typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
          InitialSolution;

        ZeroFunction zero_function(gv);
        InitialVelocity init_velocity(zero_function);
        InitialPressure init_pressure(gv);
        InitialSolution initial_solution(init_velocity,init_pressure);

        typedef BCTypeParam_PressureDrop<GV> BoundaryFunction;
        typedef PressureDropFlux<GV,RF> NeumannFlux;
        typedef ZeroVectorFunction<GV,RF,2> ExternalForce;

        // Domain parameters:
        const int tube_direction = 0; // Tube in x-axes direction
        const RF tube_length = 6.0;
        const RF tube_origin = -1.0;

        const RF boundary_pressure = lp.get<double>("boundaries.pressure");
        BoundaryFunction boundary_function(tube_length, tube_origin, tube_direction);
        NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("LU2",grid,level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
      }
#endif

#if HAVE_UG
    // UG Grid Hagen-Poiseuille test 3D
    if (params.hasSub("HU3") && params.get<bool>("HU3.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("HU3");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");
        merge_sub_trees(params.sub("boundaries"),lp,"boundaries");

        typedef Dune::UGGrid<3> GridType;
        GridType grid;

        typedef double RF;

        std::string grid_file = grid_base_dir + "/pipe.msh";
        Dune::GridFactory<GridType> factory(&grid);
        Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
        factory.createGrid();

        grid.globalRefine(level);

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
        typedef ZeroVectorFunction<GV,RF,3> ExternalForce;

        BoundaryFunction boundary_function;
        NeumannFlux neumann_flux(gv);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("HU3",grid,level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
      }
#endif

#if HAVE_UG
    // UG Grid turbulence tube test 3D
    if (params.hasSub("TU3") && params.get<bool>("TU3.enabled"))
      {
        Dune::ParameterTree& lp = params.sub("TU3");

        const int level = lp.get<int>("level");
        const int runs = lp.get("runs",global_runs);
        const bool solve = lp.get("solve",global_solve);
        const bool io = lp.get("io",global_io);

        merge_sub_trees(params.sub("physics"),lp,"physics");
        merge_sub_trees(params.sub("boundaries"),lp,"boundaries");

        typedef Dune::UGGrid<3> GridType;
        GridType grid;

        std::string grid_file = grid_base_dir + "/turbtube.msh";
        Dune::GridFactory<GridType> factory(&grid);
        Dune::GmshReader<GridType>::read(factory,grid_file,true,false);
        factory.createGrid();
        grid.globalRefine(level);

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
        typedef Dune::PDELab::PowerGridFunction<ZeroFunction,3>
          InitialVelocity;
        typedef ZeroFunction InitialPressure;
        typedef Dune::PDELab::CompositeGridFunction<InitialVelocity,InitialPressure>
          InitialSolution;

        ZeroFunction zero_function(gv);
        InitialVelocity init_velocity(zero_function);
        InitialPressure init_pressure(gv);
        InitialSolution initial_solution(init_velocity,init_pressure);


        typedef BCTypeParam_PressureDrop<GV> BoundaryFunction;
        typedef PressureDropFlux<GV,RF> NeumannFlux;
        typedef ZeroVectorFunction<GV,RF,3> ExternalForce;

        // Domain parameters:
        const int tube_direction = 2; // Tube in z-axes direction
        const RF tube_length = 5.0;
        const RF tube_origin = 0.0;

        const RF boundary_pressure = lp.get<double>("boundaries.pressure");
        BoundaryFunction boundary_function(tube_length, tube_origin, tube_direction);
        NeumannFlux neumann_flux(gv, boundary_pressure, tube_length, tube_origin, tube_direction);
        ExternalForce externalForce(gv);

        typedef Dune::PDELab::NavierStokesDefaultParameters<
          GV,
          RF,
          ExternalForce,
          BoundaryFunction,
          InitialVelocity,
          NeumannFlux
          > LOPParameters;
        LOPParameters parameters(lp.sub("physics"),
                                 externalForce,
                                 boundary_function,
                                 init_velocity,
                                 neumann_flux);

        // solve problem
        navierstokes<GV,V_FEM,P_FEM,InitialSolution,LOPParameters,q>
          (gv,name("TU3",grid,level,solve,io),parameters,vFem,pFem,initial_solution,solve,io,runs);
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
