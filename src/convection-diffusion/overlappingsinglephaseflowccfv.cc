// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Problems A-F in parallel using cell-centered finite volumes (works on nonoverlapping grids in overlapping mode)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/diffusionccfv.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"../utility/gridexamples.hh"
#include"problemA.hh"
#include"problemB.hh"
#include"problemC.hh"
#include"problemD.hh"
#include"problemE.hh"
#include"problemF.hh"

template<class GV>
void test (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim)); // works only for cubes

  // make function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::P0ParallelConstraints,VBE> GFS;
  watch.reset();
  GFS gfs(gv,fem);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // local operator
  watch.reset();
  typedef k_A<GV,RF> KType;
  KType k(gv);
  // Dune::FieldVector<double,dim> correlation_length;
  // correlation_length = 1.0/64.0;
  // KType k(gv,correlation_length,0.5,0.0,5000,-1083);
  typedef A0_A<GV,RF> A0Type;
  A0Type a0(gv);
  typedef F_A<GV,RF> FType;
  FType f(gv);
  typedef B_A<GV> BType;
  BType b(gv);
  typedef J_A<GV,RF> JType;
  JType j(gv);
  typedef G_A<GV,RF> GType;
  GType g(gv);
  typedef Dune::PDELab::DiffusionCCFV<KType,A0Type,FType,BType,JType,GType> LOP;
  LOP lop(k,a0,f,b,j,g);
  std::cout << "=== local operator setup " <<  watch.elapsed() << " s" << std::endl;

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::constraints(gfs,cc,false);

  // grid operator

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeros per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs);
  x = 0.0;
  Dune::PDELab::interpolate(g,gfs,x);

  // typedef  Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GO> LS;
  // LS ls(gfs,5000,3);

  typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<GFS,CC> LS;
  LS ls(gfs,cc,100,5,2);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x,1e-10);
  slp.apply();

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::nonconforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.pwrite("single_phase_yasp2d_CCFV","vtk","",Dune::VTK::appendedraw);
}

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

    if (argc!=4)
      {
        std::cout << "usage: " << argv[0] << " <nx> <ny> <nz>" << std::endl;
        return 0;
      }
    int nx; sscanf(argv[1],"%d",&nx);
    int ny; sscanf(argv[2],"%d",&ny);
    int nz; sscanf(argv[3],"%d",&nz);

    // 2D
    if (true)
    {
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(128));
      N[0] = nx; N[1] = ny;
      std::bitset<dim> B(false);
      int overlap=3;
      Dune::YaspGrid<dim> grid(helper.getCommunicator(),L,N,B,overlap);
      //      grid.globalRefine(6);

      // solve problem :)
      test(grid.leafGridView());
    }

    // Q1, 3d
    if (false)
    {
      // make grid
      const int dim = 3;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N;
      N[0] = nx; N[1] = ny; N[2] = nz;
      std::bitset<dim> B(false);
      int overlap=1;
      Dune::YaspGrid<dim> grid(helper.getCommunicator(),L,N,B,overlap);

      // solve problem :)
      test(grid.leafGridView());
    }

    // UG Q1 2D test
// #if HAVE_UG
//     if (false)
//     {
//       // make grid
//       UGUnitSquareQ grid(1000);
//       grid.globalRefine(6);

//       test(grid.leafGridView());
//     }
// #endif

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
