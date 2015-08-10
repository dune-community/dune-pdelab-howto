// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Problems A-F using lowest order Raviart-Thomas elements (sequential, works only with a direct solver only)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/grid/yaspgrid.hh>

#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif

#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#endif

#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/raviartthomas0.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/permeability_adapter.hh>
#include<dune/pdelab/localoperator/diffusionmixed.hh>


//===============================================================
// Choose among one of the problems A-F here:
//===============================================================
#include "parameterC.hh"
#define PARAMETERCLASS ParameterC
#define PROBLEMNAME "C"

//===============================================================
// dummy boundary condition function for the pressure component
//===============================================================

class BCTypeParam_Dummy
  : public Dune::PDELab::DirichletConstraintsParameters /*@\label{bcp:base}@*/
{
public:

  template<typename I>
  bool isDirichlet(
                   const I & intersection,   /*@\label{bcp:name}@*/
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return false;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

template<typename GV,
         typename PFEM,
         typename VFEM,
         typename PROBLEM>
void driver( const GV& gv,
             const PFEM& pfem,
             const VFEM& vfem,
             PROBLEM& problem,
             std::string filename )
{
  // types and constants
  typedef typename PFEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType Real;


  // make a grid function space
  typedef Dune::PDELab::istl::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,PFEM,Dune::PDELab::NoConstraints,VBE> P0GFS;
  P0GFS p0gfs(gv,pfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,VFEM,Dune::PDELab::RT0Constraints,VBE> RT0GFS;
  RT0GFS rt0gfs(gv,vfem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::LexicographicOrderingTag,
    RT0GFS,
    P0GFS> MGFS;
  MGFS mgfs(rt0gfs,p0gfs); // the mixed grid function space

  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> BCType;
  BCType bctype(gv,problem);

  BCTypeParam_Dummy d;

  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> GType;
  GType g(gv,problem);

  // construct a composite boundary condition type function
  typedef Dune::PDELab::CompositeConstraintsParameters<BCType,BCTypeParam_Dummy> BCT;
  BCT bct(bctype,d);

  // constraints
  typedef typename MGFS::template ConstraintsContainer<Real>::Type T;
  T t;                               // container for transformation
  Dune::PDELab::constraints(bct,mgfs,t); // fill container

  // construct a composite grid function
  typedef Dune::PDELab::ConvectionDiffusionVelocityExtensionAdapter<PROBLEM> VType;
  VType v(gv,problem);

  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rv(v);
  typedef Dune::PDELab::CompositeGridFunction<RVType,GType> UType;
  UType u(rv,g);

  // make coefficent Vectors
  typedef typename Dune::PDELab::BackendVectorSelector<MGFS,Real>::Type X;
  X x(mgfs,0.0);

  // do interpolation
  // Dune::PDELab::interpolate(u,mgfs,x);
  // Dune::PDELab::set_nonconstrained_dofs(t,0.0,x);  // clear interior

  // make grid operator
  typedef Dune::PDELab::DiffusionMixed<PROBLEM> LOP;
  LOP lop(problem,4,2);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeros per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<MGFS,MGFS,LOP,MBE,Real,Real,Real,T,T> GO;
  GO go(mgfs,t,mgfs,t,lop,mbe);

  // represent operator as a matrix
  typedef typename GO::Jacobian M;
  M m(go);
  std::cout << m.patternStatistics() << std::endl;
  m = 0.0;
  go.jacobian(x,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // set up solver
  typedef typename M::BaseT ISTLM;
#if HAVE_SUPERLU
  Dune::SuperLU<ISTLM> solver(Dune::PDELab::Backend::native(m), true);
  Dune::InverseOperatorResult stat;

  X r(mgfs,0.0);
  go.residual(x,r);
  X z(mgfs,0.0);
  solver.apply(z,r,stat);
  x -= z;
#else
#error No superLU support, please install and configure it.
#endif

  // select subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS,Dune::TypeTree::TreePath<0> > VSUB;
  VSUB vsub(mgfs);                   // velocity subspace

  typedef Dune::PDELab::GridFunctionSubSpace<MGFS,Dune::TypeTree::TreePath<1> > PSUB;
  PSUB psub(mgfs);                   // pressure subspace

  // make discrete function objects
  typedef PermeabilityAdapter<PROBLEM> PermDGF;
  PermDGF permdgf(gv,problem);
  typedef Dune::PDELab::VTKGridFunctionAdapter<PermDGF> PermVTKDGF;

  typedef Dune::PDELab::DiscreteGridFunctionPiola<VSUB,X> RT0DGF;
  RT0DGF rt0dgf(vsub,x);

  typedef Dune::PDELab::DiscreteGridFunction<PSUB,X> P0DGF;
  P0DGF p0dgf(psub,x);

  // output grid function with VTKWriter
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1); // plot result
  vtkwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P0DGF> >(p0dgf,"pressure"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<RT0DGF> >(rt0dgf,"velocity"));
  vtkwriter.addCellData(std::make_shared<PermVTKDGF>(permdgf,"logK"));
  vtkwriter.write(filename,Dune::VTK::ascii);
  std::cout << "View result using:\n paraview --data=" << filename << ".vtu \n" << std::endl;
}



int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // read command line arguments
    if (argc!=4) {
      std::cout << "usage: rt0main <dim> <geometry> <maxlevel>" << std::endl;
      std::cout << "       <dim> = 2 | 3" << std::endl;
      std::cout << "       <mesh> = cube | simplex" << std::endl;
      std::cout << "       <maxlevel> = a nonnegative integer" << std::endl;
      std::cout << std::endl;
      std::cout << "e.g.: ./rt0main 2 cube 5" << std::endl;
      std::cout << std::endl;
      return 0;
    }


    typedef double Real;
    int dim_dyn; sscanf(argv[1],"%d",&dim_dyn);
    std::string mesh(argv[2]);
    int maxlevel; sscanf(argv[3],"%d",&maxlevel);

    // 2D YASP
    if (mesh=="cube" && dim_dyn==2) {

      const int dim = 2;
      Dune::FieldVector<Real,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(maxlevel);

      // instantiate finite element maps
      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV gv = grid.leafGridView();

      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,Real,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,Real,0,Dune::GeometryType::cube> RT0FEM;
      RT0FEM rt0fem(gv);

      std::stringstream filename;
      filename << "rt0q_" << PROBLEMNAME << "_Yasp2d_level" << maxlevel;
      typedef PARAMETERCLASS<GV,Real> Problem;
      Problem problem(gv);
      driver(gv,p0fem,rt0fem,problem,filename.str());

    }

    // 3D YASP
    if (mesh=="cube" && dim_dyn==3) {

      const int dim = 3;
      Dune::FieldVector<Real,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(maxlevel);

      // instantiate finite element maps
      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV gv = grid.leafGridView();

      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,Real,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,Real,0,Dune::GeometryType::cube> RT0FEM;
      RT0FEM rt0fem(gv);

      std::stringstream filename;
      filename << "rt0q_" << PROBLEMNAME << "_Yasp3d_level" << maxlevel;
      typedef PARAMETERCLASS<GV,Real> Problem;
      Problem problem(gv);
      driver(gv,p0fem,rt0fem,problem,filename.str());

    }

#if HAVE_DUNE_ALUGRID
    // 2D ALUGRID
    if (mesh=="simplex" && dim_dyn==2) {

      const int dim = 2;
      // make grid
      typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> Grid;
      Dune::FieldVector<Grid::ctype, Grid::dimension> ll(0.0);
      Dune::FieldVector<Grid::ctype, Grid::dimension> ur(1.0);
      std::array<unsigned int, Grid::dimension> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      grid->globalRefine(maxlevel);

      // instantiate finite element maps
      typedef Grid::LeafGridView GV;
      const GV gv = grid->leafGridView();

      typedef Dune::PDELab::P0LocalFiniteElementMap<Grid::ctype, Real, dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV, Grid::ctype, Real, 0, Dune::GeometryType::simplex> RT0FEM;

      RT0FEM rt0fem(gv);

      std::stringstream filename;
      filename << "rt0q_" << PROBLEMNAME << "_Alu2d_level" << maxlevel;
      typedef PARAMETERCLASS<GV,Real> Problem;
      Problem problem(gv);
      driver(gv,p0fem,rt0fem,problem,filename.str());

    }
#endif

#if HAVE_UG
    if (mesh=="simplex" && dim_dyn==2) {
      // make grid
      const int dim = 2;
      typedef Dune::UGGrid<2> Grid;
      Dune::FieldVector<Grid::ctype, dim> ll(0.0);
      Dune::FieldVector<Grid::ctype, dim> ur(1.0);
      std::array<unsigned int, dim> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      grid->globalRefine(maxlevel);

      // instantiate finite element maps
      typedef Grid::LeafGridView GV;
      const GV gv = grid->leafGridView();

      typedef Dune::PDELab::P0LocalFiniteElementMap<Grid::ctype, Real, dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV, Grid::ctype, Real, 0, Dune::GeometryType::simplex> RT0FEM;
      RT0FEM rt0fem(gv);

      std::stringstream filename;
      filename << "rt0q_" << PROBLEMNAME << "_Ug2d_level" << maxlevel;
      typedef PARAMETERCLASS<GV,Real> Problem;
      Problem problem(gv);
      driver(gv,p0fem,rt0fem,problem,filename.str());

    }
#endif

#if HAVE_ALBERTA
    // 2D ALBERTA
    if (mesh=="simplex" && dim_dyn==2) {

      // make grid
      const int dim = 2;
      typedef Dune::AlbertaGrid<dim, dim> Grid;
      Dune::FieldVector<Grid::ctype, Grid::dimension> ll(0.0);
      Dune::FieldVector<Grid::ctype, Grid::dimension> ur(1.0);
      std::array<unsigned int, Grid::dimension> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      grid->globalRefine(maxlevel);

      // instantiate finite element maps
      typedef Grid::LeafGridView GV;
      const GV gv = grid->leafGridView();

      typedef Dune::PDELab::P0LocalFiniteElementMap<Grid::ctype, Real, dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV, Grid::ctype, Real, 0, Dune::GeometryType::simplex> RT0FEM;
      RT0FEM rt0fem(gv);

      std::stringstream filename;
      filename << "rt0q_" << PROBLEMNAME << "_Alberta2d_level" << maxlevel;
      typedef PARAMETERCLASS<GV,Real> Problem;
      Problem problem(gv);
      driver(gv,p0fem,rt0fem,problem,filename.str());

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
