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
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include"../utility/gridexamples.hh"

#define PROBLEM_A

#ifdef PROBLEM_A
#include "parameterA.hh"
#endif

#ifdef PROBLEM_B
#include "parameterB.hh"
#endif

#ifdef PROBLEM_C
#include "parameterC.hh"
#endif

#ifdef PROBLEM_D
#include "parameterD.hh"
#endif

#ifdef PROBLEM_E
#include "parameterE.hh"
#endif

#ifdef PROBLEM_F
#include "parameterF.hh"
#endif


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
           double alpha
           )
{
  // coordinate and result type
  typedef double Real;
  const int dim = GV::Grid::dimension;

  std::stringstream fullname;
  fullname << basename << "_" << method << "_w" << weights << "_k" << degree << "_dim" << dim << "_level" << level;

  // make grid function space
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE;
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
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<UDGF>(udgf,"u_h"));
      //vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<G>(g,"u"));
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

  // read command line arguments
  if (argc!=6)
    {
      std::cout << "usage: diffusion <mesh> <dim> <level> <method> <degree>" << std::endl;
      std::cout << "       <mesh> = cube | simplex" << std::endl;
      std::cout << "       <dim> = 2 | 3" << std::endl;
      std::cout << "       <level> = a nonnegative integer" << std::endl;
      std::cout << "       <method> = SIPG" << std::endl;
      std::cout << "       <degree> : polynomial degree (integer)" << std::endl;
      return 0;
    }

  std::string mesh(argv[1]);
  int dim_dyn; sscanf(argv[2],"%d",&dim_dyn);
  int level; sscanf(argv[3],"%d",&level);
  std::string method(argv[4]);
  int degree_dyn; sscanf(argv[5],"%d",&degree_dyn);

  try
    {
      if (mesh=="cube" && dim_dyn==2)
        {
          const int dim = 2;
          Dune::FieldVector<double,dim> L(1.0);
          Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
          std::bitset<dim> P(false);
          typedef Dune::YaspGrid<dim> Grid;
          Grid grid(L,N,P,0);
          typedef Grid::LeafGridView GV;

          // refine grid
          grid.globalRefine( level );
            {
              const GV& gv=grid.leafGridView();

#ifdef PROBLEM_A
              typedef ParameterA<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("A");
#endif
#ifdef PROBLEM_B
              typedef ParameterB<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("B");
#endif
#ifdef PROBLEM_C
              typedef ParameterC<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("C");
#endif
#ifdef PROBLEM_D
              typedef ParameterD<GV,double> PROBLEM;
              Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
              correlation_length = 1.0/64.0;
              PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
              std::string problemlabel("D");
#endif
#ifdef PROBLEM_E
              typedef ParameterE<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("E");
#endif
#ifdef PROBLEM_F
              typedef ParameterF<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("F");
#endif
              problemlabel.append("_CUBE");

              if (method=="SIPG") {
                if (degree_dyn==1) {
                  const int degree=1;
                  typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==2) {
                  const int degree=2;
                  typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==3) {
                  const int degree=3;
                  typedef Dune::PDELab::QkDGLocalFiniteElementMap<Grid::ctype,double,degree,dim> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
              }
            }
        }
      if (mesh=="cube" && dim_dyn==3)
        {
          const int dim = 3;
          Dune::FieldVector<double,dim> L(1.0);
          Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
          std::bitset<dim> P(false);
          typedef Dune::YaspGrid<dim> Grid;
          Grid grid(L,N,P,0);
          typedef Grid::LeafGridView GV;

          // refine grid
          grid.globalRefine( level );
            {
              const GV& gv=grid.leafGridView();

#ifdef PROBLEM_A
              typedef ParameterA<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("A");
#endif
#ifdef PROBLEM_B
              typedef ParameterB<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("B");
#endif
#ifdef PROBLEM_C
              typedef ParameterC<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("C");
#endif
#ifdef PROBLEM_D
              typedef ParameterD<GV,double> PROBLEM;
              Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
              correlation_length = 1.0/64.0;
              PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
              std::string problemlabel("D");
#endif
#ifdef PROBLEM_E
              typedef ParameterE<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("E");
#endif
#ifdef PROBLEM_F
              typedef ParameterF<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("F");
#endif
              problemlabel.append("_CUBE");

              if (method=="SIPG") {
                if (degree_dyn==1) {
                  const int degree=1;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);

                }
                if (degree_dyn==2) {
                  const int degree=2;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==3) {
                  const int degree=3;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==4) {
                  const int degree=4;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
              }
            }
        }

#if HAVE_ALUGRID
      if (mesh=="simplex" && dim_dyn==2)
        {
          const int dim = 2;

          // make grid
          ALUUnitCube<dim> unitcube;
          typedef ALUUnitCube<dim>::GridType Grid;
          typedef Grid::LeafGridView GV;

          unitcube.grid().globalRefine( level );
            {
              const GV& gv=unitcube.grid().leafView();

#ifdef PROBLEM_A
              typedef ParameterA<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("A");
#endif
#ifdef PROBLEM_B
              typedef ParameterB<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("B");
#endif
#ifdef PROBLEM_C
              typedef ParameterC<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("C");
#endif
#ifdef PROBLEM_D
              typedef ParameterD<GV,double> PROBLEM;
              Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
              correlation_length = 1.0/64.0;
              PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
              std::string problemlabel("D");
#endif
#ifdef PROBLEM_E
              typedef ParameterE<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("E");
#endif
#ifdef PROBLEM_F
              typedef ParameterF<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("F");
#endif

              problemlabel.append("_SIMPLEX");

              if (method=="SIPG") {
                if (degree_dyn==1) {
                  const int degree=1;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==2) {
                  const int degree=2;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==3) {
                  const int degree=3;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
              }
            }
        }
      if (mesh=="simplex" && dim_dyn==3)
        {
          const int dim = 3;

          // make grid
          ALUUnitCube<dim> unitcube;
          typedef ALUUnitCube<dim>::GridType Grid;
          typedef Grid::LeafGridView GV;


          unitcube.grid().globalRefine( level );
            {
              const GV& gv=unitcube.grid().leafView();

#ifdef PROBLEM_A
              typedef ParameterA<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("A");
#endif
#ifdef PROBLEM_B
              typedef ParameterB<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("B");
#endif
#ifdef PROBLEM_C
              typedef ParameterC<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("C");
#endif
#ifdef PROBLEM_D
              typedef ParameterD<GV,double> PROBLEM;
              Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
              correlation_length = 1.0/64.0;
              PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
              std::string problemlabel("D");
#endif
#ifdef PROBLEM_E
              typedef ParameterE<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("E");
#endif
#ifdef PROBLEM_F
              typedef ParameterF<GV,double> PROBLEM;
              PROBLEM problem;
              std::string problemlabel("F");
#endif

              problemlabel.append("_SIMPLEX");

              if (method=="SIPG") {
                if (degree_dyn==1) {
                  const int degree=1;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==2) {
                  const int degree=2;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
                if (degree_dyn==3) {
                  const int degree=3;
                  typedef Dune::PDELab::OPBLocalFiniteElementMap<Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
                  FEMDG femdg;
                  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
                  runDG<GV,FEMDG,PROBLEM,degree,blocksize>(gv,femdg,problem,problemlabel,level,"SIPG","ON",2.0);
                }
              }
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
