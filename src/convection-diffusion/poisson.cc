// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Poisson problem on various grids (sequential).
    HANGING_NODES_REFINEMENT is macro used to switch on hanging nodes tests.
    It is set in "Makefile.am" to generate the executable 'poisson_HN'.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/float_cmp.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/hangingnode.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/backend/seqistlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include "../utility/gridexamples.hh"

/*
  HANGING_NODES_REFINEMENT is macro used to switch on hanging nodes tests.
  It is set in "Makefile.am" to generate the executable 'poisson_HN'.
*/

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

// function for defining the source term
template<typename GV, typename RF>
class F
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F<GV,RF> > BaseT;

  F (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    if (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375)
      y = 50.0;
    else
      y = 0.0;
    y=0;
  }
};



// constraints parameter class for selecting boundary condition type
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters /*@\label{bcp:base}@*/
{
public:

  template<typename I>
  bool isDirichlet(
				   const I & intersection,   /*@\label{bcp:name}@*/
				   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
				   ) const
  {

    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );

    if( xg[1]<1E-6 || xg[1]>1.0-1E-6 )
      return false; // Neumann b.c.
    else if( xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6 )
      return false; // Neumann b.c.
    else
      return true;  // Dirichlet b.c. on all other boundaries
  }

  template<typename I>
  bool isNeumann(
                 const I & intersection,   /*@\label{bcp:name}@*/
                 const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                 ) const
  {
    return !isDirichlet(intersection,coord);
  }

};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G<GV,RF> > BaseT;

  G (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    typename Traits::DomainType center;
    for (int i=0; i<GV::dimension; i++) center[i] = 0.5;
    center -= x;
	y = exp(-center.two_norm2());
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J<GV,RF> > BaseT;

  J (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
							  typename Traits::RangeType& y) const
  {
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      {
        y = 0;
        return;
      }
    if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
      {
        y = -5.0;
        return;
      }
  }
};







//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename BCTYPE, typename CON>
void poisson_driver(const GV& gv,
                    const FEM& fem,
                    std::string filename,
                    const BCTYPE& bctype,    // boundary condition type
                    bool hanging_nodes,
                    int q,  // quadrature order
                    const CON& con = CON())
{
  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make grid function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);
  gfs.name("poisson solution");

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();

  Dune::PDELab::constraints(bctype,gfs,cg);

  // make grid operator
  typedef F<GV,R> FType;
  FType f(gv);
  typedef J<GV,R> JType;
  JType j(gv);
  typedef Dune::PDELab::Poisson<FType,BCTypeParam,JType> LOP;
  LOP lop(f,bctype,j,q);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(45); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,R,R,R,C,C> GO;
  GO go(gfs,cg,gfs,cg,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,R> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);

  // Choose ISTL Solver Backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
  LS ls(5000,2);

  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x0,1e-12);
  slp.setHangingNodeModifications(hanging_nodes);
  slp.apply();

  // make discrete function object
  Dune::SubsamplingVTKWriter<GV> vtkwriter( gv, 1 );
  //Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x0);
  vtkwriter.write(filename,Dune::VTK::ascii);
}







#ifdef HANGING_NODES_REFINEMENT

template<typename Grid>
void doSomeRandomRefinement( Grid & grid ){

  // Do some random refinement. The result is a grid that may
  // contain multiple hanging nodes per edge!

  typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
    Iterator;

  for(int i=0; i<4;++i){
    Iterator it = grid.template leafbegin<0,Dune::All_Partition>();
    Iterator eit = grid.template leafend<0,Dune::All_Partition>();

    for(;it!=eit;++it){
      if((double)rand()/(double)RAND_MAX > 0.6)
        grid.mark(1,*(it));
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
  }

}

#endif









//===============================================================
// Main program with grid setup
//===============================================================
// The domain is always the unit square in 2D or the unit cube in 3D
//
// Overview:
//
// Testcase 1.) ALUGrid 2D triangular cells (hanging nodes refinement) - Pk elements
// Testcase 2.) ALUGrid 3D cubical cells (hanging nodes refinement) - Q1 elements
// Testcase 3.) ALUGrid 3D tetrahedral cells (hanging nodes refinement) - P1 elements
//
// Testcase 4.) YaspGrid 2D rectangular cells (uniform refinement) - Q1 elements
// Testcase 5.) YaspGrid 2D rectangular cells (uniform refinement) - Q2 elements
// Testcase 6.) YaspGrid 3D rectangular cells (uniform refinement) - Q1 elements
//
// Testcase  7.) UG 2D triangular cells (hanging nodes refinement) - P1 elements
// Testcase  8.) UG 2D rectangular cells (hanging nodes refinement) - Q1 elements
// Testcase  9.) UG 3D cubical cells (hanging nodes refinement) - Q1 elements
// Testcase 10.) UG 3D tetrahedral cells (hanging nodes refinement) - P1 elements
//
// Testcase 11.) Alberta Grid 2D triangular cells (uniform refinement) - Pk elements
//
// Not supported by the Grid: ALUGrid 2D rectangular cells

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_ALUGRID
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 1.) ALUGrid 2D triangular cells (hanging nodes refinement) - P1 elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Testcase 1.) ALUGrid 2D triangular cells (uniform refinement) - P1 elements"
        << std::endl;
#endif

      // make grid
      typedef ALUUnitSquare Grid;
      Grid grid;
      grid.globalRefine(4);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=1; // polynomial order of the FEM
      const int q=2*k; // integration order for the quadrature rule
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      // We need the boundary function for the hanging nodes
      // constraints engine as we have to distinguish between hanging
      // nodes on dirichlet and on neumann boundaries
      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      // This is the type of the local constraints assembler that has
      // to be adapted for different local basis spaces and grid types
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::SimplexGridP1Assembler ConstraintsAssembler;

      // The type of the constraints engine
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;

      // Get constraints engine. We set adaptToIsolateHangingNodes =
      // true and therefore the constructor refines the grid until
      // there are less than two hanging nodes per edge.
      Constraints constraints(grid,true,bctype);

      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_ALU_Pk_2d_randomly_refined",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_ALU_Pk_2d",bctype,false,q);
#endif
    }
#endif


#if HAVE_ALUGRID
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 2.) ALUGrid 3D cubical cells (hanging nodes refinement) - Q1 elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Alternative Testcase 2.) ALUGrid 3D cubical cells (uniform refinement) - Q1 elements"
        << std::endl;
#endif

      // make grid
      typedef ALUCubeUnitSquare Grid;
      Grid grid;
      grid.globalRefine(2);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=1;
      const int q=2*k;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::CubeGridQ1Assembler ConstraintsAssembler;
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;

      Constraints constraints(grid,true,bctype);

      // solve problem
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_ALU_Q1_3d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_ALU_Q1_3d",bctype,false,q);
#endif

    }
#endif


#if HAVE_ALUGRID
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 3.) ALUGrid 3D tetrahedral cells (hanging nodes refinement) - Pk elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Alternative Testcase 3.) ALUGrid 3D tetrahedral cells (uniform refinement) - Pk elements"
        << std::endl;
#endif
      // make grid
      typedef ALUUnitCube<3> UnitCube;
      UnitCube unitcube;
      typedef ALUUnitCube<3>::GridType Grid;
      Grid &grid = unitcube.grid();
      grid.globalRefine(1);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map

      typedef UnitCube::GridType::ctype DF;
      typedef double R;
      const int k=1;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      BCTypeParam bctype;


#ifdef HANGING_NODES_REFINEMENT
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::SimplexGridP1Assembler ConstraintsAssembler;
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;
      Constraints constraints(grid,true,bctype);
      // solve problem
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_ALU_Pk_3d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>( gv,fem,"poisson_ALU_Pk_3d", bctype, false, q );
#endif


    }
#endif


    //return 0;

    {
      std::cout
        << std::endl << std::endl
        << "Testcase 4.) YaspGrid 2D rectangular cells (uniform refinement) - Q1 elements"
        << std::endl;

      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      std::bitset<2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(6);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      const int k=1;
      const int q=2*k;
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      BCTypeParam bctype;
      // solve problem
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_yasp_Q1_2d",bctype,false,q);
    }


    {
      std::cout
        << std::endl << std::endl
        << "Testcase 5.) YaspGrid 2D rectangular cells (uniform refinement) - Q2 elements"
        << std::endl;

      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      std::bitset<2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      const int k=2;
      const int q=2*k;
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
      BCTypeParam bctype;

      // solve problem
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_yasp_Q2_2d",bctype,false,q);
    }

    {
      std::cout
        << std::endl << std::endl
        << "Testcase 6.) YaspGrid 3D rectangular cells (uniform refinement) - Q1 elements"
        << std::endl;

      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      std::bitset<3> B(false);
      Dune::YaspGrid<3> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      const int k=2;
      const int q=2*k;
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
      BCTypeParam bctype;

      // solve problem
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_yasp_Q1_3d",bctype,false,q);
    }

#if HAVE_UG
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 7.) UG 2D triangular cells (hanging nodes refinement) - P1 elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Alternative Testcase 7.) UG 2D triangular cells (uniform refinement) - P1 elements"
        << std::endl;
#endif

      // make grid
      typedef UGUnitSquare Grid;
      Grid grid;
      grid.setRefinementType( Grid::LOCAL );
      grid.setClosureType( Grid::NONE );  // This is needed to get hanging nodes refinement! Otherwise you would get triangles.
      grid.globalRefine(4);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=1; //k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      // We need the boundary function for the hanging nodes
      // constraints engine as we have to distinguish between hanging
      // nodes on dirichlet and on neumann boundaries
      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::SimplexGridP1Assembler ConstraintsAssembler;
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;
      Constraints constraints(grid,true,bctype);
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_UG_Pk_2d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_UG_Pk_2d",bctype,false,q);
#endif

    }
#endif // HAVE_UG





#if HAVE_UG
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 8.) UG 2D rectangular cells (hanging nodes refinement) - Q1 elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Alternative Testcase 8.) UG 2D rectangular cells (uniform refinement) - Q1 elements"
        << std::endl;
#endif
      // make grid (unitcube made of cubes)
      typedef UGUnitSquareQ Grid;
      Grid grid;

      grid.setRefinementType( Grid::LOCAL );
      grid.setClosureType( Grid::NONE );  // This is needed to get hanging nodes refinement! Otherwise you would get triangles.
      grid.globalRefine(4);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      const int k=1;
      const int q=2*k;
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::CubeGridQ1Assembler ConstraintsAssembler;
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;
      Constraints constraints(grid,true,bctype);
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_UG_Q1_2d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_UG_Q1_2d",bctype,false,q);
#endif

    }
#endif // HAVE_UG






#ifdef HAVE_UG
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 9.) UG 3D cubical cells (hanging nodes refinement) - Q1 elements"
        << std::endl;
#endif
      // get grid and do a single global refine
      typedef UGUnitCube<3,1>::GridType Grid;
      UGUnitCube<3,1> ugunitcube;
      Grid & grid = ugunitcube.grid();
      grid.setRefinementType(Grid::LOCAL);
      grid.setClosureType(Grid::NONE);
      grid.globalRefine(1);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;

      // make finite element map
      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      const GV& gv=grid.leafGridView();

      const int k=1;
      const int q=2*k;

      typedef GV::Grid::ctype DF;
      typedef double R;

      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::CubeGridQ1Assembler ConstraintsAssembler;

      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;

      Constraints constraints(grid,true,bctype);
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_UG_Q1_3d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#endif

    }
#endif// HAVE_UG






#ifdef HAVE_UG
    {
#ifdef HANGING_NODES_REFINEMENT
      std::cout
        << std::endl << std::endl
        << "Testcase 10.) UG 3D tetrahedral cells (hanging nodes refinement) - P1 elements"
        << std::endl;
#else
      std::cout
        << std::endl << std::endl
        << "Alternative Testcase 10.) UG 3D tetrahedral cells (uniform refinement) - P1 elements"
        << std::endl;
#endif

      // UG Grid made of tetrahedrons - test Pk3D with hanging nodes!
      typedef UGUnitCube<3,2>::GridType Grid;
      UGUnitCube<3,2> ugunitcube;
      Grid & grid = ugunitcube.grid();
      grid.setRefinementType(Grid::LOCAL);
      grid.setClosureType(Grid::NONE);
      grid.globalRefine(1);

#ifdef HANGING_NODES_REFINEMENT
      doSomeRandomRefinement<Grid>( grid );
#endif
      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=1;     // polynomial degree
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      BCTypeParam bctype;

#ifdef HANGING_NODES_REFINEMENT
      typedef Dune::PDELab::HangingNodesConstraintsAssemblers::SimplexGridP1Assembler ConstraintsAssembler;
      typedef Dune::PDELab::HangingNodesDirichletConstraints
        <GV::Grid,ConstraintsAssembler,BCTypeParam> Constraints;
      Constraints constraints(grid,true,bctype);
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,
                                                       fem,
                                                       "poisson_UG_Pk_3d_hangingNodes",
                                                       bctype,
                                                       true,
                                                       q,
                                                       constraints
                                                       );
#else
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_UG_Pk_3d",bctype,false,q);
#endif

    }
#endif // HAVE_UG








#if HAVE_ALBERTA
    {
      std::cout
        << std::endl << std::endl
        << "Testcase 11.) Alberta 2D triangular cells (uniform refinement) - Pk elements"
        << std::endl;

      // make grid
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
      BCTypeParam bctype;

      // solve problem
      typedef Dune::PDELab::ConformingDirichletConstraints Constraints;
      poisson_driver<GV,FEM,BCTypeParam,Constraints>(gv,fem,"poisson_Alberta_Pk_2d",bctype,false,q);
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
