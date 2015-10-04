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
#include<dune/pdelab/localoperator/convectiondiffusionccfv.hh>
#include<dune/pdelab/localoperator/darcy_CCFV.hh>
#include<dune/pdelab/localoperator/permeability_adapter.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"parameterA.hh"

template<typename GV>
void test (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;

  typedef ParameterA<GV,RF> PROBLEM;
  PROBLEM problem(gv);

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

  typedef Dune::PDELab::ConvectionDiffusionCCFV<PROBLEM> LOP;
  LOP lop(problem);

  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);

  std::cout << "=== local operator setup " <<  watch.elapsed() << " s" << std::endl;

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::constraints(gfs,cc,false);

  // grid operator
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(2*dim+1); // Maximal number of nonzeros per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs);
  x = 0.0;
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  // make vector consistent
  {
      Dune::PDELab::CopyDataHandle<GFS,V> dh(gfs,x);
      if(gfs.gridView().comm().size() > 1)
          gfs.gridView().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

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
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF> >(dgf,"solution"));

  typedef DarcyVelocityFromHeadCCFV<PROBLEM,DGF> DarcyDGF;
  DarcyDGF darcydgf(problem,dgf);
  typedef Dune::PDELab::VTKGridFunctionAdapter<DarcyDGF> DarcyVTKDGF;
  vtkwriter.addVertexData(std::make_shared<DarcyVTKDGF>(darcydgf,"darcyvelocity"));

  typedef PermeabilityAdapter<PROBLEM> PermDGF;
  PermDGF permdgf(gv,problem);
  typedef Dune::PDELab::VTKGridFunctionAdapter<PermDGF> PermVTKDGF;
  vtkwriter.addCellData(std::make_shared<PermVTKDGF>(permdgf,"logK"));

  std::stringstream filename;
  filename << "single_phase_yasp" << dim << "d_CCFV";
  vtkwriter.pwrite(filename.str(),"vtk","",Dune::VTK::appendedraw);
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
      Dune::YaspGrid<dim> grid(L,N,B,overlap);
      //      grid.globalRefine(6);

      // solve problem :)
      test(grid.leafGridView());
    }

    // Q1, 3d
    if (true)
    {
      // make grid
      const int dim = 3;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N;
      N[0] = nx; N[1] = ny; N[2] = nz;
      std::bitset<dim> B(false);
      int overlap=1;
      Dune::YaspGrid<dim> grid(L,N,B,overlap);

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
