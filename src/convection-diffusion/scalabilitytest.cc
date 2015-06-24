// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief High-level test with Poisson equation
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
#include<dune/common/typetraits.hh>
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

#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/localoperator/convectiondiffusionccfv.hh>
#include<dune/pdelab/localoperator/darcy_CCFV.hh>
#include<dune/pdelab/localoperator/permeability_adapter.hh>

#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include "parameter_factory.hh"


const bool graphics = true;

template<typename GV,typename PROBLEM>
void test_ccfv (const GV& gv,PROBLEM& problem)
{
  typedef typename GV::Grid::ctype DF;
  typedef typename PROBLEM::RangeFieldType RF;
  const int dim = GV::dimension;

  std::stringstream fullname;
  fullname << "scalabilitytest_CCFV_dim" << dim;

  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim)); // works only for cubes

  // make function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::P0ParallelConstraints,VBE> GFS;
  watch.reset();
  GFS gfs(gv,fem);

  // local operator
  watch.reset();

  typedef Dune::PDELab::ConvectionDiffusionCCFV<PROBLEM> LOP;
  LOP lop(problem);

  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::constraints(g,gfs,cc,false);

  // grid operator
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
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
  LS ls(gfs,cc,10,5,1);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x,1e-6);
  slp.apply();

  // make discrete function object
  if( graphics ){

    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> UDGF;
    UDGF udgf(gfs,x);
    Dune::VTKWriter<GV> vtkwriter(gv);
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF> >(udgf,"pressure_h"));

    typedef DarcyVelocityFromHeadCCFV<PROBLEM,UDGF> DarcyDGF;
    DarcyDGF darcydgf(problem,udgf);
    typedef Dune::PDELab::VTKGridFunctionAdapter<DarcyDGF> DarcyVTKDGF;
    vtkwriter.addVertexData(std::make_shared<DarcyVTKDGF>(darcydgf,"velocity_h"));

    typedef PermeabilityAdapter<PROBLEM> PermDGF;
    PermDGF permdgf(gv,problem);
    typedef Dune::PDELab::VTKGridFunctionAdapter<PermDGF> PermVTKDGF;
    vtkwriter.addCellData(std::make_shared<PermVTKDGF>(permdgf,"logK"));

    vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);

  }

}



//! solve problem with DG method
template<class GV, class FEM, class PROBLEM, int degree, int blocksize>
void runDG ( const GV& gv,
             const FEM& fem,
             PROBLEM& problem,
             std::string basename,
             int level,
             std::string method,
             std::string weights,
             double alpha )
{
  // coordinate and result type
  typedef double Real;
  const int dim = GV::Grid::dimension;

  std::stringstream fullname;
  fullname << "scalabilitytest_";

  fullname << basename << "_" << method << "_w" << weights << "_k" << degree << "_dim" << dim << "_level" << level;

  // make grid function space
  typedef Dune::PDELab::P0ParallelConstraints CON;
  const Dune::PDELab::ISTLParameters::Blocking blocking
    = Dune::PDELab::ISTLParameters::static_blocking;
  typedef Dune::PDELab::ISTLVectorBackend<blocking,blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // make local operator
  Dune::PDELab::ConvectionDiffusionDGMethod::Type m;
  if (method=="SIPG") m = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
  if (method=="NIPG") m = Dune::PDELab::ConvectionDiffusionDGMethod::NIPG;
  Dune::PDELab::ConvectionDiffusionDGWeights::Type w;
  if (weights=="ON") w = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;
  if (weights=="OFF") w = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOff;
  typedef Dune::PDELab::ConvectionDiffusionDG<PROBLEM,FEM> LOP;
  LOP lop(problem,m,w,alpha);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(g,gfs,cc,false);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make a vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);

  // make linear solver and solve problem
  int verbose=1;
  if (gv.comm().rank()!=0) verbose=0;
  if (method=="SIPG")
    {
      typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<GFS,CC> LS;
      LS ls(gfs,cc,10,5,verbose);
      // typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
      // LS ls(10000,1);
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      SLP slp(go,ls,u,1e-6);
      slp.apply();
    }
  else
    {
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
      LS ls(gfs,cc,10,5,verbose);
      // typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
      // LS ls(10000,1);
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      SLP slp(go,ls,u,1e-6);
      slp.apply();
    }

  if( graphics ){
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> UDGF;
    UDGF udgf(gfs,u);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,std::max(0,degree-1));
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF> >(udgf,"u_h"));
    vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);
  }

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

  if (argc!=9 && argc!=6) {
    if(helper.rank()==0) {
      std::cout << "usage option 1: " << argv[0] << " <problem> <order> <nx> <ny> <nz>" << std::endl;
      std::cout << "usage option 2: " << argv[0] << " <problem> <order> <nx> <ny> <nz> <px> <py> <pz>" << std::endl;
      std::cout << "       <problem> = A | C | E " << std::endl;
    }
    return 0;
  }
  char choice; sscanf(argv[1],"%c",&choice);
  int degree_dyn; sscanf(argv[2],"%d",&degree_dyn);
  int nx; sscanf(argv[3],"%d",&nx);
  int ny; sscanf(argv[4],"%d",&ny);
  int nz; sscanf(argv[5],"%d",&nz);

  int px=0; int py=0; int pz=0;
  if (argc==9){
    sscanf(argv[6],"%d",&px);
    sscanf(argv[7],"%d",&py);
    sscanf(argv[8],"%d",&pz);
  }

  try
    {
      typedef double Real;

      const int dim = 3;
      Dune::FieldVector<Real,dim> L(1.0);
      Dune::array<int,dim> N;
      N[0] = nx; N[1] = ny; N[2] = nz;
      std::bitset<dim> B(false);
      int overlap=1;

      typedef Dune::YaspFixedSizePartitioner<dim> YP;
      YP* yp = (YP*) Dune::YaspGrid<dim>::defaultLoadbalancer();
      if( px*py*pz==0 ){
        // If px,py,pz were not specified choose the default load balancer
        if( helper.rank() == 0 )
          std::cout << "Using default partitioning of YASP." << std::endl;
      }

      else if( px*py*pz != helper.size() ){
        // If px*py*pz is not equal to the available number of processors
        // wrong input, stop and output warning!
        if( helper.rank()==0 )
          std::cerr << "Wrong input: px*py*pz != np" << std::endl;
        exit(1);
      }

      else {
        std::array<int,dim> yasppartitions;
        yasppartitions[0] = px;
        yasppartitions[1] = py;
        yasppartitions[2] = pz;
        yp = new YP(yasppartitions);
      }

      Dune::YaspGrid<dim> grid(L,N,B,overlap,helper.getCommunicator(),yp);

      typedef Dune::YaspGrid<dim> Grid;
      typedef Grid::LeafGridView GV;

      const GV& gv=grid.leafGridView();

      typedef ParameterBase<GV,Real> PROBLEM;
      typedef ParameterFactory<PROBLEM,char> ProblemFactory;
      ProblemFactory::template registerAll<GV,Real,char>(gv);
      PROBLEM* pProblem = ProblemFactory::getInstance().createParameter( gv, choice );

      if (degree_dyn==0) {
        //test_old_ccfv(gv);
        test_ccfv(gv,*pProblem);
      }
      if (degree_dyn==1) {
        const int degree=1;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,"CUBE",0,"SIPG","ON",2.0);
      }
      if (degree_dyn==2) {
        const int degree=2;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,"CUBE",0,"SIPG","ON",2.0);
      }
      if (degree_dyn==3) {
        const int degree=3;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,"CUBE",0,"SIPG","ON",2.0);
      }
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
