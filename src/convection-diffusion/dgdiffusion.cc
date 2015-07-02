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
#include<dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#endif

#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include "parameter_factory.hh"


const bool graphics = true;


/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                   typename T1::Traits::RangeFieldType,
                   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceSquaredAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor
  DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename T1::Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename T2::Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1.two_norm2();
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }

private:
  const T1& t1;
  const T2& t2;
};

//! solve problem with DG method
template<class GV, class FEM, class PROBLEM, int degree, int blocksize>
void runDG(
           const GV& gv,
           const FEM& fem,
           PROBLEM& problem,
           std::string basename,
           int level,
           std::string method,
           std::string weights,
           typename PROBLEM::RangeFieldType alpha
           )
{
  // coordinate and result type
  typedef typename PROBLEM::RangeFieldType Real;
  const int dim = GV::Grid::dimension;

  std::stringstream fullname;
  fullname << "vtk/dgdiffusion_" << basename << "_" << method << "_w" << weights << "_k" << degree << "_dim" << dim << "_level" << level;

  // make grid function space
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,blocksize> VBE;
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
  MBE mbe(9); // number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make a vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);

  // make linear solver and solve problem
  int ls_verbosity = 2;
  if (method=="SIPG")
    {
      typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
      LS ls(10000,ls_verbosity);

      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      SLP slp(go,ls,u,1e-12);
      slp.apply();
    }
  else
    {
      //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GFS> LS;
      //LS ls(2,10000,ls_verbosity);
      typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
      LS ls(10000,ls_verbosity);
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      SLP slp(go,ls,u,1e-12);
      slp.apply();
    }

  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> UDGF;
  UDGF udgf(gfs,u);

  // compute L2 error if analytical solution is available
  //typedef DifferenceSquaredAdapter<G,UDGF> DifferenceSquared;
  //DifferenceSquared differencesquared(g,udgf);
  //typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  //Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,12);
  //std::cout << fullname.str()
  //          << " N=" << std::setw(11) << gfs.globalSize()
  //          << " L2ERROR=" << std::setw(11) << std::setprecision(3) << std::scientific << std::uppercase << sqrt(l2errorsquared[0]) << std::endl;

  // write vtk file
  if (graphics)
    {
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,degree-1);
      vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<UDGF>>(udgf,"u_h"));
      //vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<G> >(g,"u"));
      vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);
      std::cout << "Run: \n paraview --data=" << fullname.str() << ".vtu \n" << std::endl;
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

  // read command line arguments
  if (argc!=7)
    {
      std::cout << "usage: diffusion <problem> <dim> <geometry> <level> <method> <degree>" << std::endl;
      std::cout << "       <problem>  = A | B | ... | F" << std::endl;
      std::cout << "       <dim>      = 2 | 3" << std::endl;
      std::cout << "       <geometry> = cube | simplex" << std::endl;
      std::cout << "       <level>    = global refinement level > 0" << std::endl;
      std::cout << "       <method>   = SIPG | NIPG" << std::endl;
      std::cout << "       <degree>   = polynomial degree (integer)" << std::endl;
      std::cout << std::endl;
      std::cout << "e.g.: ./dgdiffusion C 2 cube 3 SIPG 1" << std::endl;
      std::cout << std::endl;
      return 0;
    }

  typedef double Real;

  char choice; sscanf(argv[1],"%c",&choice);
  int dim_dyn; sscanf(argv[2],"%d",&dim_dyn);
  std::string mesh(argv[3]);
  int level; sscanf(argv[4],"%d",&level);
  std::string method(argv[5]);
  int degree_dyn; sscanf(argv[6],"%d",&degree_dyn);

  try {

    // 2D YASP:
    if (mesh=="cube" && dim_dyn==2) {
      const int dim = 2;
      Dune::FieldVector<Real,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> P(false);
      typedef Dune::YaspGrid<dim> Grid;
      Grid grid(L,N,P,0);
      typedef Grid::LeafGridView GV;
      typedef ParameterBase<GV,Real> PROBLEM;
      typedef ParameterFactory<PROBLEM,char> ProblemFactory;

      // refine grid
      grid.globalRefine( level );

      const GV& gv=grid.leafGridView();
      ProblemFactory::template registerAll<GV,Real,char>(gv);
      PROBLEM* pProblem = ProblemFactory::getInstance().createParameter( gv, choice );
      std::string problemlabel("X");
      problemlabel[0] = choice;
      problemlabel.append("_CUBE");

      if (degree_dyn==1) {
        const int degree=1;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==2) {
        const int degree=2;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==3) {
        const int degree=3;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,Real,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }

    }


    // 3D YASP:
    if (mesh=="cube" && dim_dyn==3) {
      const int dim = 3;
      Dune::FieldVector<Real,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> P(false);
      typedef Dune::YaspGrid<dim> Grid;
      Grid grid(L,N,P,0);
      typedef Grid::LeafGridView GV;
      typedef ParameterBase<GV,Real> PROBLEM;
      typedef ParameterFactory<PROBLEM,char> ProblemFactory;

      // refine grid
      grid.globalRefine( level );

      const GV& gv=grid.leafGridView();
      ProblemFactory::template registerAll<GV,Real,char>(gv);
      PROBLEM* pProblem = ProblemFactory::getInstance().createParameter( gv, choice );
      std::string problemlabel("X");
      problemlabel[0] = choice;
      problemlabel.append("_CUBE");

      if (degree_dyn==1) {
        const int degree=1;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::cube> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);

      }
      if (degree_dyn==2) {
        const int degree=2;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::cube> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==3) {
        const int degree=3;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::cube> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==4) {
        const int degree=4;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::cube> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }

    }


#if HAVE_DUNE_ALUGRID
    if (mesh=="simplex" && dim_dyn==2) {
      // make grid'
      const int dim = 2;
      typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> Grid;
      Dune::FieldVector<Grid::ctype, Grid::dimension> ll(0.0);
      Dune::FieldVector<Grid::ctype, Grid::dimension> ur(1.0);
      std::array<unsigned int, Grid::dimension> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      typedef Grid::LeafGridView GV;
      typedef ParameterBase<GV,Real> PROBLEM;
      typedef ParameterFactory<PROBLEM,char> ProblemFactory;

      grid->globalRefine( level );

      GV gv = grid->leafView();
      ProblemFactory::template registerAll<GV,Real,char>(gv);
      PROBLEM* pProblem = ProblemFactory::getInstance().createParameter( gv, choice );
      std::string problemlabel("X");
      problemlabel[0] = choice;
      problemlabel.append("_SIMPLEX");

      if (degree_dyn==1) {
        const int degree=1;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==2) {
        const int degree=2;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==3) {
        const int degree=3;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }

    }

    // 3D ALUGRID
    if (mesh=="simplex" && dim_dyn==3) {
      // make grid
      const int dim = 3;
      typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> Grid;
      Dune::FieldVector<Grid::ctype, Grid::dimension> ll(0.0);
      Dune::FieldVector<Grid::ctype, Grid::dimension> ur(1.0);
      std::array<unsigned int, Grid::dimension> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      typedef Grid::LeafGridView GV;
      typedef ParameterBase<GV,Real> PROBLEM;
      typedef ParameterFactory<PROBLEM,char> ProblemFactory;

      grid->globalRefine( level );

      GV gv = grid->leafView();
      ProblemFactory::template registerAll<GV,Real,char>(gv);
      PROBLEM* pProblem = ProblemFactory::getInstance().createParameter( gv, choice );
      std::string problemlabel("X");
      problemlabel[0] = choice;
      problemlabel.append("_SIMPLEX");

      if (degree_dyn==1) {
        const int degree=1;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==2) {
        const int degree=2;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }
      if (degree_dyn==3) {
        const int degree=3;
        typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,Real,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,*pProblem,problemlabel,level,method,"ON",2.0);
      }

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
