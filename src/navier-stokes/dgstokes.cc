// -*- tab-width: 2; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Stokes with DG method (stationary case).
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
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
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
#include <dune/pdelab/localoperator/dgnavierstokes.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include "navierstokes_initial.hh"

#define USE_SUPER_LU
#define MAKE_VTK_OUTPUT

//===============================================================
// Problem setup and solution
//===============================================================

// generate a composite P2/P1 function and output it
template<typename GV, typename RF, int vOrder, int pOrder>
void stokes (const GV& gv, const Dune::ParameterTree& configuration, std::string filename)
{
  // <<<1>>> constants and types
  using ES = Dune::PDELab::AllEntitySet<GV>;
  ES es(gv);
  typedef typename ES::Grid::ctype DF;
  static const unsigned int dim = ES::dimension;
  Dune::Timer watch;
  std::cout << "=== Initialize" << std::endl;

  // <<<2>>> Make grid function space
  watch.reset();
  typedef Dune::PDELab::MonomLocalFiniteElementMap<DF,RF,dim,vOrder> vFEM;
  typedef Dune::PDELab::MonomLocalFiniteElementMap<DF,RF,dim,pOrder> pFEM;

  vFEM vFem(Dune::GeometryType(Dune::GeometryType::cube,dim));
  pFEM pFem(Dune::GeometryType(Dune::GeometryType::cube,dim));
  // DOFs per cell
  static const unsigned int vBlockSize = Dune::MonomImp::Size<dim,vOrder>::val;
  static const unsigned int pBlockSize = Dune::MonomImp::Size<dim,pOrder>::val;

  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::none,vBlockSize> VVectorBackend;
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::none,pBlockSize> PVectorBackend;
  typedef Dune::PDELab::istl::VectorBackend<> VelocityVectorBackend;

#if 1
  // this creates a flat backend (i.e. blocksize == 1)
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::none> VectorBackend;
#else
  // this creates a backend with static blocks matching the size of the LFS
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed> VectorBackend;
#endif
  // velocity
  Dune::dinfo << "--- v^dim" << std::endl;
  typedef Dune::PDELab::EntityBlockedOrderingTag VelocityOrderingTag;
  typedef Dune::PDELab::VectorGridFunctionSpace<
    ES,vFEM,dim,
    VelocityVectorBackend,
    VVectorBackend,
    Dune::PDELab::NoConstraints,
    VelocityOrderingTag
    > velocityGFS;
  velocityGFS velocityGfs(es,vFem);
  velocityGfs.name("v");
  // p
  Dune::dinfo << "--- p" << std::endl;
  typedef Dune::PDELab::GridFunctionSpace<ES,pFEM,
                                          Dune::PDELab::NoConstraints, PVectorBackend> pGFS;
  pGFS pGfs(es,pFem);
  pGfs.name("p");
  // GFS
  Dune::dinfo << "--- v^dim,p" << std::endl;
  typedef Dune::PDELab::EntityBlockedOrderingTag StokesOrderingTag;
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VectorBackend, StokesOrderingTag,
    velocityGFS, pGFS> GFS;
  GFS gfs(velocityGfs, pGfs);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // <<<3>>> Make coefficient Vector and initialize it from a function
  using V = Dune::PDELab::Backend::Vector<GFS,RF>;
  V x(gfs);

  typedef BCTypeParamGlobalDirichlet BType;
  BType b;
  typedef ZeroVectorFunction<ES,RF,dim> FType;
  FType f(es);
  typedef HagenPoiseuilleVelocityBox<ES,RF,dim> VType;
  VType v(es);
  typedef ZeroScalarFunction<ES,RF> PType;
  PType p(es);

  // <<<4>>> Make grid Function operator
  watch.reset();
  typedef Dune::PDELab::DefaultInteriorPenalty<RF> PenaltyTerm;

  typedef Dune::PDELab::DGNavierStokesParameters<ES,RF,FType,BType,VType,PType,false,false,PenaltyTerm>
    LocalDGOperatorParameters;
  LocalDGOperatorParameters lop_params(configuration.sub("parameters"),f,b,v,p);
  typedef Dune::PDELab::DGNavierStokes<LocalDGOperatorParameters> LocalDGOperator;
  const int superintegration_order = 0;
  LocalDGOperator lop(lop_params,superintegration_order);

  typedef Dune::PDELab::EmptyTransformation C;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(75); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator
    <GFS,GFS,LocalDGOperator,MBE,RF,RF,RF,C,C> GOS;
  GOS gos(gfs,gfs,lop,mbe);

  std::cout << "=== grid operator space setup " <<  watch.elapsed() << " s" << std::endl;

  typedef typename GOS::Jacobian M;
  watch.reset();
  M m(gos);
  std::cout << m.patternStatistics() << std::endl;
  std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
  m = 0.0;
  watch.reset();
  gos.jacobian(x,m);
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;

  // std::ofstream matrix("Matrix");
  // Dune::printmatrix(matrix, m.base(), "M", "r", 6, 3);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  watch.reset();
  x = 0.0;
  gos.residual(x,r);
  std::cout << "=== residual evaluation " <<  watch.elapsed() << " s" << std::endl;

  bool verbose = true;

  using Dune::PDELab::Backend::native;
  typedef Dune::PDELab::Backend::Native<M> ISTLM;
  typedef Dune::PDELab::Backend::Native<V> ISTLV;
#ifdef USE_SUPER_LU // use lu decomposition as solver
#if HAVE_SUPERLU
  // make ISTL solver
  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
  Dune::SuperLU<ISTLM> solver(native(m), verbose?1:0);
  Dune::InverseOperatorResult stat;
#else
#error No superLU support, please install and configure it.
#endif
#else // Use iterative solver
  // make ISTL solver
  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(native(m));
  Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(native(m),1.0);
  Dune::BiCGSTABSolver<ISTLV> solver(opa,ilu0,1E-10,20000, verbose?2:1);
  Dune::InverseOperatorResult stat;
#endif

  // solve the jacobian system
  r *= -1.0; // need -residual
  x = r;
  solver.apply(native(x),native(r),stat);

  //    #ifdef MAKE_VTK_OUTPUT
  // output grid function with SubsamplingVTKWriter
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x);
  vtkwriter.write(filename,Dune::VTK::appendedraw);
  //    #endif
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    int x=2;
    int y=2;
    int z=2;
    if (argc > 1)
      x = atoi(argv[1]);
    if (argc > 2)
      y = atoi(argv[2]);
    if (argc > 3)
      z = atoi(argv[3]);

    Dune::ParameterTree configuration;
    const std::string config_filename("dgstokes.ini");
    std::cout << "Reading ini-file \""<< config_filename
              << "\"" << std::endl;

    Dune::ParameterTreeParser::readINITree(config_filename, configuration);

    // YaspGrid P2/P1 2D test
    if(1){
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N;
      N[0] = x; N[1] = y;
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);

      // get view
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV gv=grid.leafGridView();

      // solve problem
      Dune::dinfo.push(false);
      stokes<GV,double,2,1>(gv,configuration,"dgstokes-2D-2-1");
    }


    // YaspGrid P2/P1 3D test
#if 1
    {
      // make grid
      const int dim = 3;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N;
      N[0] = x; N[1] = y; N[2] = z;
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);

      // get view
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV gv=grid.leafGridView();

      // solve problem
      stokes<GV,double,2,1>(gv,configuration,"dgstokes-3D-2-1");
    }
#endif

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
