// -*- tab-width: 2; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Navier-Stokes with DG method and vector valued velocity function space
    (stationary case).
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#endif
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/dgnavierstokesvelvecfem.hh>
#include <dune/pdelab/finiteelementmap/brezzidouglasmarinifem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/newton/newton.hh>

#include "navierstokes_initial.hh"

//===============================================================
// The driver for all examples
//===============================================================

template<typename GV, typename RF, typename vFEM, typename pFEM>
void navierstokesvecfem(
                        const GV& gv,
                        std::string filename,
                        vFEM vFem,
                        pFEM pFem)
{
  // constants and types
  typedef typename GV::Grid::ctype DF;
  const unsigned int dim = GV::dimension;

  // vector backend for velocity and pressure
  typedef Dune::PDELab::ISTLVectorBackend<> VelocityVectorBackend;
  typedef Dune::PDELab::ISTLVectorBackend<> PVectorBackend;

  // this creates a flat backend (i.e. blocksize == 1)
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::none> VectorBackend;

  // velocity grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,vFEM,Dune::PDELab::NoConstraints,VelocityVectorBackend> velocityGFS;
  velocityGFS velocityGfs(gv,vFem);
  velocityGfs.name("v");

  // pressure grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,pFEM,Dune::PDELab::NoConstraints,PVectorBackend> pGFS;
  pGFS pGfs(gv,pFem);
  pGfs.name("p");

  // composite grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VectorBackend,
    Dune::PDELab::EntityBlockedOrderingTag,
    velocityGFS,
    pGFS>
    GFS;
  GFS gfs(velocityGfs,pGfs);

  // boundary conditions and parameters
  typedef BCTypeParamGlobalDirichlet BType;
  BType b;
  typedef ZeroVectorFunction<GV,RF,dim> FType;
  FType f(gv);
  typedef HagenPoiseuilleVelocityBox<GV,RF,dim> VType;
  VType v(gv);
  typedef ZeroScalarFunction<GV,RF> PType;
  PType p(gv);

  const RF mu = 1.0;
  const RF rho = 1.0;
  std::string method = "-1 20.0 1";
  typedef typename Dune::PDELab::DefaultInteriorPenalty<RF> PenaltyTerm;
  PenaltyTerm ip_term(method,mu);

  typedef Dune::PDELab::DGNavierStokesParameters<GV,RF,FType,BType,VType,PType,true,false,PenaltyTerm>
    LocalDGOperatorParameters;
  LocalDGOperatorParameters lop_params(method,mu,rho,f,b,v,p,ip_term);

  // local operators
  typedef Dune::PDELab::DGNavierStokesVelVecFEM<LocalDGOperatorParameters> LocalDGOperator;
  LocalDGOperator lop(lop_params,1);

  // grid operators
  typedef Dune::PDELab::EmptyTransformation C;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(75); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  typedef Dune::PDELab::GridOperator<
    GFS,GFS,LocalDGOperator,MBE,RF,RF,RF,C,C> GOS;
  GOS gos(gfs,gfs,lop,mbe);

  typedef typename GOS::Traits::Domain V;
  typedef typename GOS::Traits::Jacobian M;

  // random initial guess
  V x(gfs);
  for(typename V::size_type i=0; i<x.N(); i++)
    Dune::PDELab::Backend::native(x)[i] = (static_cast<RF>(random()))/RAND_MAX;

  // linear solver
  typedef Dune::PDELab::ISTLBackend_SEQ_GMRES_ILU0 LinearSolver;
  LinearSolver ls;

  // solve problem
  Dune::PDELab::Newton<GOS,LinearSolver,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(30);
  newton.apply();

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write(filename,Dune::VTK::appendedraw);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  //Maybe initialize Mpi
  Dune::MPIHelper::instance(argc, argv);

  std::string example_switch;
  if(argc != 2) {
    std::cout << std::endl << "PDELab Navier-Stokes DG VectorFEM examples" << std::endl
              << "--------------------------------------------" << std::endl << std::endl
              << "Call with command line parameter to execute examples:" << std::endl << std::endl
              << "            Setup                        | Parameter " << std::endl
              << "-----------------------------------------------------" << std::endl
              << "Hagen-Poiseuille 2D - YaspGrid   - BDM1/Q0 : HY2 " << std::endl
              << "Hagen-Poiseuille 2D - AluSimplex - BDM1/P0 : HA2 " << std::endl
              << std::endl << std::endl
              << "You might also want to take a look at the configuration file \"navierstokes_initial.hh\"."
              << std::endl << std::endl;
    exit(1);
  }
  else
    example_switch = argv[1];

  try {

    // YaspGrid Hagen-Poiseuille test
    if(example_switch.find("HY2") != std::string::npos) {
      typedef double RF;
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(40));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);

      // get view
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      const int vOrder = 1;
      const int pOrder = vOrder - 1;

      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap<GV,DF,RF,vOrder> vFEM;
      vFEM vFem(gv);
      typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,pOrder,dim> pFEM;
      pFEM pFem;

      // solve problem
      navierstokesvecfem<GV,RF,vFEM,pFEM>
        (gv,"hagenpoiseuille_yasp_BDM1Q0_2d",vFem,pFem);
    }

    //==================================================================//
    // NOTE
    // -----
    // At the moment BDM1 on triangles seems to be broken.
    // Therefore we have the switch below deactivated.
    //==================================================================//

#if 0
#if HAVE_DUNE_ALUGRID
    // ALU Grid Hagen-Poiseuille test
    if(example_switch.find("HA2") != std::string::npos) {
      typedef double RF;
      // make grid
      const int dim = 2;
      typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming> Grid;
      Dune::FieldVector<Grid::ctype, Grid::dimension> ll(0.0);
      Dune::FieldVector<Grid::ctype, Grid::dimension> ur(1.0);
      std::array<unsigned int, Grid::dimension> elements;
      std::fill(elements.begin(), elements.end(), 20);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);

      // get view
      typedef Grid::LeafGridView GV;
      // TODO Why not const reference?
      const GV& gv = grid->leafGridView();

      // make finite element map
      const int vOrder = 1;
      const int pOrder = vOrder - 1;

      typedef Grid::ctype DF;
      // This doesn't work yet !!!
      typedef Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap<GV,DF,RF,vOrder> vFEM;

      // typedef Dune::PDELab::BrezziDouglasMariniLocalFiniteElementMap<GV,DF,RF,vOrder,Dune::GeometryType::simplex> vFEM;
      vFEM vFem(gv);
      typedef Dune::PDELab::MonomLocalFiniteElementMap<DF,RF,dim,pOrder> pFEM;
      pFEM pFem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

      // solve problem
      navierstokesvecfem<GV,RF,vFEM,pFEM>
        (gv,"hagenpoiseuille_alu_BDM1P0_2d",vFem,pFem);
    }
#endif
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
