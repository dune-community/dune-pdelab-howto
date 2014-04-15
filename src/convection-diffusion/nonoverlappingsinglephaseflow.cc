// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Problems A-F in parallel on non-overlapping grids using conforming linear finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

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

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template< typename PROBLEM, typename GV, typename FEM>
void driver(PROBLEM& problem, const GV& gv, const FEM& fem,
            std::string filename, int intorder=1 )
{
  // constants and types and global variables
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints<GV> CON;
  typedef Dune::PDELab::GridFunctionSpace
    <GV,FEM,CON,Dune::PDELab::ISTLVectorBackend<> > GFS;
  CON con(gv);
  GFS gfs(gv,fem,con);
  con.compute_ghosts(gfs);

  // make constraints map and initialize it from a function and ghost
  typedef typename GFS::template ConstraintsContainer<R>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> bctype(gv,problem);

  // make grid operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<PROBLEM,FEM> LOP;
  LOP lop(problem);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator
      <GFS,GFS,LOP,
       MBE,
       R,R,R,CC,CC,true> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::constraints(bctype,gfs,cc,false);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS;
  LS ls (go,5000,3,2);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x,1e-12);
  slp.apply();

  // output solution and data decomposition with VTKWriter
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.write(filename,Dune::VTK::appendedraw);
}

//===============================================================
// Main program with grid setup
//===============================================================

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

    // Q1, 2d
    if (true)
    {
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(32));
      std::bitset<dim> B(false);
      int overlap=0;
      Dune::YaspGrid<dim> grid(helper.getCommunicator(),L,N,B,overlap);
      //grid.globalRefine(4);
      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

        typedef Dune::YaspGrid<dim>::ctype DF;
        const int degree = 1;
        typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,degree> FEM;
        FEM fem(gv);

        driver(problem,gv,fem,"yasp2d_Q1",2*degree);
    }

    // Q1, 3d
    if (false)
    {
      // make grid
      const int dim = 3;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(8));
      std::bitset<dim> B(false);
      int overlap=0;
      Dune::YaspGrid<dim> grid(helper.getCommunicator(),L,N,B,overlap);
      //grid.globalRefine(3);

      typedef Dune::YaspGrid<dim>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

        typedef Dune::YaspGrid<dim>::ctype DF;
        const int degree = 1;
        typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,degree> FEM;
        FEM fem(gv);

        driver(problem,gv,fem,"yasp3d_Q1",2*degree);
    }

#if HAVE_UG
    if (false)
    {
      typedef Dune::UGGrid<3> GridType;
      GridType grid;

      // read gmsh file
      Dune::GridFactory<GridType> factory(&grid);
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;
      Dune::GmshReader<GridType>::read(factory,"grids/cube1045.msh",true,true);
      factory.createGrid();
      grid.loadBalance();

      std::cout << " after load balance /" << helper.rank() << "/ " << grid.size(0) << std::endl;
      //grid.globalRefine(1);
      // std::cout << " after refinement /" << helper.rank() << "/ " << grid.size(0) << std::endl;

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafGridView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=1;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      //driver(problem,gv,fem,"UG3d_P1",q);
    }
#endif

#if HAVE_ALUGRID
    // ALU Pk 3D test
    if (false)
    {
      typedef Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming> GridType;

      // read gmsh file
      Dune::GridFactory<GridType> factory(helper.getCommunicator());
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;
      if (helper.rank()==0)
        {
          Dune::GmshReader<GridType>::read(factory,"grids/cube1045.msh",true,false);
        }
      GridType *grid=factory.createGrid();
      grid->loadBalance();
      std::cout << " after load balance /" << helper.rank() << "/ " << grid->size(0) << std::endl;
      //grid->globalRefine(1);
      std::cout << " after refinement /" << helper.rank() << "/ " << grid->size(0) << std::endl;

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid->leafGridView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=1;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      driver(problem,gv,fem,"ALU3d_P1_PlastkDoeddel",q);
    }
	// ALU Q1 3D test
	if (false)
	{
	  // make grid
      typedef Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming> GridType;
	  Dune::GridFactory<GridType> factory(helper.getCommunicator());

	  if (helper.rank()==0)
		{
		  Dune::FieldVector< double, 3 > pos;
		  pos[0] = 0;  pos[1] = 0;	pos[2] = 0;	   factory.insertVertex(pos);
		  pos[0] = 1;  pos[1] = 0;	pos[2] = 0;	   factory.insertVertex(pos);
		  pos[0] = 0;  pos[1] = 1;	pos[2] = 0;	   factory.insertVertex(pos);
		  pos[0] = 1;  pos[1] = 1;	pos[2] = 0;	   factory.insertVertex(pos);
		  pos[0] = 0;  pos[1] = 0;	pos[2] = 1;	   factory.insertVertex(pos);
		  pos[0] = 1;  pos[1] = 0;	pos[2] = 1;	   factory.insertVertex(pos);
		  pos[0] = 0;  pos[1] = 1;	pos[2] = 1;	   factory.insertVertex(pos);
		  pos[0] = 1;  pos[1] = 1;	pos[2] = 1;	   factory.insertVertex(pos);

		  const Dune::GeometryType type( Dune::GeometryType::cube, 3 );
		  std::vector< unsigned int > cornerIDs( 8 );
		  for( int i = 0; i < 8; ++i )
			cornerIDs[ i ] = i;
		  factory.insertElement( type, cornerIDs );
		}

	  GridType *grid=factory.createGrid();

	  std::cout << " after reading /" << helper.rank() << "/ " << grid->size(0) << std::endl;
	  grid->loadBalance();
	  std::cout << " after load balance /" << helper.rank() << "/ " << grid->size(0) << std::endl;

      grid->globalRefine(4);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid->leafGridView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

	  // make finite element map
	  typedef GridType::ctype DF;
      const int degree = 1;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,degree> FEM;
      FEM fem(gv);

      driver(problem,gv,fem,"ALU3d_Q1",2*degree);
	}
#endif

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
