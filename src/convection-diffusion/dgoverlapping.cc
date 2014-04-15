// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief test overlapping Schwarz solver for discontinuous Galerkin (solves Problem A-F)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/parallel/mpihelper.hh>
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
#include<dune/istl/superlu.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/paamg/smoother.hh>
#include<dune/istl/paamg/construction.hh>
#include<dune/istl/paamg/parameters.hh>
#include<dune/istl/overlappingschwarz.hh>
#include<dune/istl/ilusubdomainsolver.hh>

#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/diffusiondg.hh>

#include"../utility/gridexamples.hh"

// Select Problem
#include"problemA.hh"  // exp(-norm(x,y))
#include"problemB.hh"  // Like problem A but corners have small parts with Neumann boundary
#include"problemC.hh"  // Constant flow with checkerboard changing permeability, kind of ground water problem
#include"problemD.hh"  // Constant flow with randomly changig permeability
#include"problemE.hh"  // Constant flow with constant permeability with is 1E-6 in x direction
#include"problemF.hh"  // Constant flow with constant permeability with is 1E-6 in x direction

// Define most changed values
#define QUADRATURE_RULE_ORDER 6
#define MONOM_BASIS_ORDER 2
#define BLOCK_SIZE 6
#define GRID_REFINE 1
#define DG_METHOD 2  // OBB: 0, NIPG: 1, SIPG: 2
#define MAKE_VTK_OUTPUT
//#define CALCULATE_L2_ERROR
//#define CALCULATE_ABSOLUTE_ERROR
#define USE_SUPER_LU
//#define REFINE_STEPWISE

template<class M>
struct AllStrong
  : public Dune::Amg::Parameters
{
  typedef M matrix_type;

  template<class M1>
  void init(const M1* m)
  {}
  template<class M1>
  void initRow(const M1& m, int i)
  {}
  template<class M1>
  void examine(const M1& m)
  {}

  template<class G, class E, class T>
  void examine(G& g, const E& e, const T& t)
  {
    e.properties().setDepends();
    e.properties().setInfluences();
  }
  bool isIsolated()
  {
    return false;
  }
  AllStrong(const Dune::Amg::Parameters& parms)
    : Dune::Amg::Parameters(parms)
  {}
  AllStrong()
    : Dune::Amg::Parameters()
  {}

};

template<class M>
struct AllStrongCriterion
  : public Dune::Amg::AggregationCriterion<AllStrong<M> >
{
public:
  AllStrongCriterion(const Dune::Amg::Parameters& parms)
    : Dune::Amg::AggregationCriterion<AllStrong<M> >(parms)
  {}
  AllStrongCriterion()
  {}
  typedef M Matrix;
};

// Solve the given problem with OBB, NIPG or SIPG
template<class GV, class FEM>
void solve_dg (const GV& gv, const FEM& fem, std::string filename, const bool verbose)
{
  typedef double RF;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,BLOCK_SIZE> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
                                          Dune::PDELab::NoConstraints,VBE> GFS;
  watch.reset();
  GFS gfs(gv,fem);
  if (verbose)
    {
      std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;
    }

  // make coefficient Vector and initialize it from a function
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
  V x0(gfs);
  x0 = 0.0;
  typedef K_A<GV,RF> KType;
  KType k(gv);
  typedef F_A<GV,RF> FType;
  FType f(gv);

  typedef BCTypeParam_A BCType;
  BCType bctype;

  typedef G_A<GV,RF> GType;
  GType g(gv);
  typedef J_A<GV,RF> JType;
  JType j(gv);
  Dune::PDELab::interpolate(g,gfs,x0);

  // make grid function operator
  typedef Dune::PDELab::DiffusionDG<KType,FType,BCType,GType,JType> LOP;
  LOP la(k,f,bctype,g,j,DG_METHOD);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF> GOS;
  GOS gos(gfs,gfs,la,mbe);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<RF>::Type M;
  watch.reset();
  M m(gos);
  if (verbose)
    {
      std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
    }
  m = 0.0;
  watch.reset();
  gos.jacobian(x0,m);
  if (verbose)
    {
      std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
    }
  //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  watch.reset();
  x0 = 0.0;
  gos.residual(x0,r);
  //std::cout << "Residuenvektor" << std::endl << r << std::endl;
  if (verbose)
    {
      std::cout << "=== residual evaluation " <<  watch.elapsed() << " s" << std::endl;
    }

  std::cout<<"Matrix size is "<<m.M()<<"x"<<m.N()<<" with block entries of size "<<BLOCK_SIZE<<"x"<<BLOCK_SIZE<<std::endl;

#ifdef USE_SUPER_LU // use lu decomposition as solver
#if HAVE_SUPERLU
  // make ISTL solver
  typedef typename V::BaseT ISTLV;
  typedef typename M::BaseT ISTLM;
  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(Dune::PDELab::istl::raw(m));
  Dune::SuperLU<ISTLM> solver(Dune::PDELab::istl::raw(m), verbose?1:0);
  Dune::InverseOperatorResult stat;
#else
#error No superLU support, please install and configure it.
#endif
#else // Use iterative solver
  // make ISTL solver


  typedef typename M::BaseT ISTLM;
  typedef typename V::BaseT ISTLV;
  typedef  Dune::SeqOverlappingSchwarz<ISTLM,ISTLV,Dune::AdditiveSchwarzMode,
#if defined SUPERLU_SD && HAVE_SUPERLU
    Dune::SuperLU<ISTLM>
#else
    Dune::ILU0SubdomainSolver<ISTLM,ISTLV,ISTLV>
#endif
    >
    DGSmoother;
  /*#if defined SUPERLU_SD && HAVE_SUPERLU
    , Dune::SuperLU<ISTLM>
    #else
    #ifdef ILUN_SD
    , Dune::ILUNSubdomainSolver<ISTLM,V,V>
    #endif
    #endif
    > DGSmoother;*/

  typedef typename Dune::Amg::SmootherTraits<DGSmoother>::Arguments DGSmootherArgs;
  typename Dune::Amg::ConstructionTraits<DGSmoother>::Arguments cargs;
  DGSmootherArgs dgSmootherArgs;
  dgSmootherArgs.relaxationFactor = .8;
  dgSmootherArgs.onthefly = true; // compute decomposition on the fly
  dgSmootherArgs.overlap=
    Dune::Amg::SeqOverlappingSchwarzSmootherArgs<typename ISTLM::field_type>::none;

  if(dgSmootherArgs.overlap!=Dune::Amg::SeqOverlappingSchwarzSmootherArgs<typename ISTLM::field_type>::pairwise){

    std::cout<<"Aggregating"<<std::endl;
    typedef typename Dune::Amg::MatrixGraph<ISTLM> MatrixGraph;
    typedef typename Dune::Amg::PropertiesGraph<MatrixGraph,Dune::Amg::VertexProperties,
                                                Dune::Amg::EdgeProperties,Dune::IdentityMap,Dune::IdentityMap> PropertiesGraph;
    MatrixGraph mg(m);
    PropertiesGraph pg(mg,Dune::IdentityMap(),Dune::IdentityMap());
    Dune::Amg::AggregatesMap<typename MatrixGraph::VertexDescriptor> map(pg.maxVertex()+1);
    //Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<ISTLM,
    //Dune::Amg::RowSum> > dgcriterion;
    Dune::Amg::CoarsenCriterion<AllStrongCriterion<ISTLM> > dgcriterion;
    dgcriterion.setDefaultValuesIsotropic(2,5);
    std::cout<<" dg criterion "<<dgcriterion<<std::endl;
    int noAggregates, isoAggregates, oneAggregates, skippedAggregates;

    Dune::tie(noAggregates, isoAggregates, oneAggregates,skippedAggregates) =
      map.buildAggregates(m, pg, dgcriterion, false);
    std::cout<<"no aggregates="<<noAggregates<<" iso="<<isoAggregates<<" one="<<oneAggregates<<" skipped="<<skippedAggregates<<std::endl;
    // misuse coarsener to renumber aggregates
    Dune::Amg::IndicesCoarsener<Dune::Amg::SequentialInformation,int> renumberer;
    typedef std::vector<bool>::iterator Iterator;
    typedef Dune::IteratorPropertyMap<Iterator, Dune::IdentityMap> VisitedMap;
    std::vector<bool> excluded(m.N(), false);
    VisitedMap vm(excluded.begin(), Dune::IdentityMap());
    Dune::Amg::SequentialInformation info;
    renumberer.coarsen(Dune::Amg::SequentialInformation(), pg, vm, map, info);
    cargs.setArgs(dgSmootherArgs);
    cargs.setMatrix(m, map);
    map.free();
  }else{
    cargs.setArgs(dgSmootherArgs);
    cargs.setMatrix(m);
  }

  DGSmoother* smoother = Dune::Amg::ConstructionTraits<DGSmoother>::construct(cargs);

  Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m);
  Dune::BiCGSTABSolver<ISTLV> solver(opa,*smoother,1E-10,10,2);//verbose?1:0);
  Dune::InverseOperatorResult stat;
#endif

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  solver.apply(x,r,stat);
  x += x0;
  //std::cout << x << std::endl;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

#ifdef MAKE_VTK_OUTPUT
  // output grid function with SubsamplingVTKWriter
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"u"));
  vtkwriter.write(filename,Dune::VTK::ascii);
#endif

#ifdef ANALYTIC_SOLUTION_PROVIDED
  AnalyticSolution<GV,RF> analyticSolution(gv);
#else
  std::cout << "No analytic solution provided, can't calculate L2 error" << std::endl;
#warning No analytic solution provided, will not calculate L2 Error
#endif

#ifdef CALCULATE_L2_ERROR
#ifdef ANALYTIC_SOLUTION_PROVIDED
  double error_l2 = l2error(analyticSolution,gfs,x,QUADRATURE_RULE_ORDER);
  std::cout.precision(8);
  std::cout << "L2 error: "
            << std::setw(8) << gv.size(0) << " elements "
            << std::scientific << error_l2 << std::endl;
#else
  std::cout << "No analytic solution provided, can't calculate L2 error" << std::endl;
#warning No analytic solution provided, will not calculate L2 Error
#endif
#endif

#ifdef CALCULATE_ABSOLUTE_ERROR
#ifdef ANALYTIC_SOLUTION_PROVIDED
  Dune::FieldVector<DF,dim> evaluate_point(0.8);
  double error_absolute = absolute_error(analyticSolution,gfs,x,evaluate_point,verbose);
  std::cout.precision(8);
  std::cout << "Absolute error: "
            << std::setw(8) << gv.size(0) << " elements "
            << std::scientific << error_absolute << std::endl;
#endif
#endif
}

int main(int argc, char** argv)
{
  //Maybe initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is a sequential program." << std::endl;
  else
    {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

  std::string problem="C";

  try
    {
      // 2D
      if (true)
        {
          // make grid
          const int dim = 2;
          Dune::FieldVector<double,dim> L(1.0);  // L[0]=2.0; L[1]=1.0;
          Dune::array<int,dim> N(Dune::fill_array<int,dim>(64));
          std::bitset<dim> B(false);
          Dune::YaspGrid<dim> grid(L,N,B,0);
#ifdef REFINE_STEPWISE
          for (int i = 0; i <= GRID_REFINE; ++i)
            {
              solve_dg(grid.leafGridView(), false);
              grid.globalRefine(1);
            }
#else
          grid.globalRefine(GRID_REFINE);

          // instantiate finite element maps
          typedef Dune::PDELab::MonomLocalFiniteElementMap<double,double,2,MONOM_BASIS_ORDER> FEM;
          FEM fem(Dune::GeometryType(Dune::GeometryType::cube,2)); // works only for cubes

          // solve problem :)
          solve_dg(grid.leafGridView(),fem,"DG_Yasp_2d",true);
#endif
        }

#if HAVE_ALBERTA
      if (true)
        {
          typedef AlbertaUnitSquare GridType;
          GridType grid;
          grid.globalRefine(10);

          // get view
          typedef GridType::LeafGridView GV;
          const GV& gv=grid.leafGridView();

          // instantiate finite element maps
          typedef Dune::PDELab::MonomLocalFiniteElementMap<double,double,2,MONOM_BASIS_ORDER> FEM;
          FEM fem(Dune::GeometryType(Dune::GeometryType::simplex,2)); // works only for cubes

          solve_dg(gv,fem,"DG_Alberta_2d",true);

        }
#endif


    }
  catch (Dune::Exception &e)
    {
      std::cerr << "Dune reported error: " << e << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << "Unknown exception thrown!" << std::endl;
      return 1;
    }
}
