// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Solve Problems A-F in parallel using cell-centered finite volumes (works on nonoverlapping grids in overlapping mode)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p0constraints.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/diffusionccfv.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
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
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::P0ParallelConstraints,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS; 
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
  Dune::PDELab::constraints(g,gfs,cc,false);

  // grid operator
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,VBE::MatrixBackend,RF,RF,RF,
    CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs);
  x = 0.0;
  Dune::PDELab::interpolate(g,gfs,x);

  // typedef  Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GO> LS;
  // LS ls(gfs,5000,3);
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(gfs,cc,5000,5,1);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,x,ls,1e-8);
  slp.apply();

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.pwrite("single_phase_yasp2d_CCFV","vtk","",Dune::VTKOptions::binaryappended);
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
    // if (false)
    // {
    //   // make grid
    //   Dune::FieldVector<double,2> L(1.0);
    //   Dune::FieldVector<int,2> N(128);
    //   Dune::FieldVector<bool,2> B(false);
    //   int overlap=4;
    //   Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
    //   //      grid.globalRefine(6);
      
    //   // solve problem :)
    //   test(grid.leafView());
    // }

    // Q1, 3d
    if (true)
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N;
      N[0] = nx; N[1] = ny; N[2] = nz;
      Dune::FieldVector<bool,3> B(false);
      int overlap=1;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);

      // solve problem :)
      test(grid.leafView());
    }

    // UG Q1 2D test
// #if HAVE_UG
//     if (false)
//     {
//       // make grid 
//       UGUnitSquareQ grid(1000);
//       grid.globalRefine(6);

//       test(grid.leafView());
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
