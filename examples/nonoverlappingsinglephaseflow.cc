// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
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

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/rannacher_turek2dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/diffusion.hh>

#include"gridexamples.hh"
#include"problemA.hh"
#include"problemB.hh"
#include"problemC.hh"
#include"problemD.hh"
#include"problemE.hh"
#include"problemF.hh"

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template<typename BType, typename GType, typename KType, typename A0Type, typename FType, typename JType,
         typename GV, typename FEM, typename FEM0> 
void driver (const BType& b, const GType& g, 
             const KType& k, const A0Type& a0, const FType& f, const JType& j,
             const GV& gv, const FEM& fem, const FEM0& fem0, std::string filename, int intorder=1)
{
  // constants and types and global variables
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::NonoverlappingConformingDirichletConstraints,
    Dune::PDELab::ISTLVectorBackend<1>,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  Dune::PDELab::NonoverlappingConformingDirichletConstraints cn;
  GFS gfs(gv,fem,cn);
  cn.compute_ghosts(gfs);

  // make constraints map and initialize it from a function and ghost
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(b,gfs,cg);
  std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size() << std::endl;

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<R>::Type V; 
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // make grid function operator
  typedef Dune::PDELab::Diffusion<KType,A0Type,FType,BType,JType> LOP; 
  LOP lop(k,a0,f,b,j,intorder);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1>,true> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  watch.reset();
  M m(gos); m = 0.0;
  std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
  watch.reset();
  gos.jacobian(x,m);
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // solve the jacobian system
  V r(gfs,0.0);
  gos.residual(x,r);
  V z(gfs,0.0);

  // set up parallel solver
  typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;
  PHELPER phelper(gfs);
  typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,V> POP;
  POP pop(gfs,m);
  typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
  PSP psp(gfs,phelper);
  typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,V> PRICH;
  PRICH prich(gfs,phelper);
  int verbose=0;
  if (gv.comm().rank()==0) verbose=1;
  Dune::CGSolver<V> solver(pop,psp,prich,1E-6,250000,verbose);
  Dune::InverseOperatorResult stat;  
  solver.apply(z,r,stat);
  x -= z;

  // output solution and data decomposition with VTKWriter
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,
    Dune::PDELab::NoConstraints,
    Dune::PDELab::ISTLVectorBackend<1>,
    Dune::PDELab::SimpleGridFunctionStaticSize
    > GFS0; 
  GFS0 gfs0(gv,fem0);
  typedef typename GFS0::template VectorContainer<R>::Type V0;
  V0 partition(gfs0,0.0);
  Dune::PDELab::PartitionDataHandle<GFS0,V0> pdh(gfs0,partition);
  if (gfs.gridview().comm().size()>1)
    gfs0.gridview().communicate(pdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  typedef Dune::PDELab::DiscreteGridFunction<GFS0,V0> DGF0;
  DGF0 pdgf(gfs0,partition);

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
  //Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(pdgf,"decomposition"));
  vtkwriter.write(filename,Dune::VTKOptions::binaryappended);
}

template<typename GV, typename FEM, typename FEM0> 
void dispatcher (std::string problem, const GV& gv, const FEM& fem, const FEM0& fem0, 
                 std::string gridname, int intorder=1)
{
  std::string A("A"), B("B"), C("C"), D("D"), E("E"), F("F");
  std::string filename(""), underscore("_");
  filename = problem+underscore+gridname;

  typedef double RF;

  if (problem==A) 
    {
      driver(B_A<GV>(gv), G_A<GV,RF>(gv),K_A<GV,RF>(gv),
             A0_A<GV,RF>(gv),F_A<GV,RF>(gv),J_A<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
  if (problem==B) 
    {
      driver(B_B<GV>(gv), G_B<GV,RF>(gv),K_B<GV,RF>(gv),
             A0_B<GV,RF>(gv),F_B<GV,RF>(gv),J_B<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
  if (problem==C) 
    {
      driver(B_C<GV>(gv), G_C<GV,RF>(gv),K_C<GV,RF>(gv),
             A0_C<GV,RF>(gv),F_C<GV,RF>(gv),J_C<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
  if (problem==D) 
    {
      Dune::FieldVector<RF,GV::Grid::dimension> correlation_length;
      correlation_length = 1.0/64.0;
      driver(B_D<GV>(gv), G_D<GV,RF>(gv),
             K_D<GV,RF>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,RF>(gv),F_D<GV,RF>(gv),J_D<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
  if (problem==E) 
    {
      driver(B_E<GV>(gv), G_E<GV,RF>(gv),K_E<GV,RF>(gv),
             A0_E<GV,RF>(gv),F_E<GV,RF>(gv),J_E<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
  if (problem==F) 
    {
      driver(B_F<GV>(gv), G_F<GV,RF>(gv), K_F<GV,RF>(gv),
             A0_F<GV,RF>(gv), F_F<GV,RF>(gv), J_F<GV,RF>(gv),
             gv,fem,fem0,filename,intorder);
    }
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
    
    std::string problem="C";

#if HAVE_MPI
    // Q1, 2d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(256);
      Dune::FieldVector<bool,2> B(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
      //grid.globalRefine(4);
      
      typedef Dune::YaspGrid<2>::ctype DF;
      typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,2> FEM0;
      FEM0 fem0(Dune::GeometryType::cube);

      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();

      dispatcher(problem,gv,fem,fem0,"yasp2d_Q1",2);
    }
#endif

#if HAVE_MPI
    // Q1, 3d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(8);
      Dune::FieldVector<bool,3> B(false);
      int overlap=0;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
      //grid.globalRefine(3);

      typedef Dune::YaspGrid<3>::ctype DF;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,3> FEM0;
      FEM0 fem0(Dune::GeometryType::cube);

      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafView();

      dispatcher(problem,gv,fem,fem0,"yasp3d_Q1",2);
    }
#endif

#if HAVE_MPI
    // Rannacher Turek, 2d
    if (true)
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
      grid.globalRefine(5);
      
      typedef Dune::YaspGrid<2>::ctype DF;
      typedef Dune::PDELab::RannacherTurek2DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,double,2> FEM0;
      FEM0 fem0(Dune::GeometryType::cube);

      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();

      dispatcher(problem,gv,fem,fem0,"yasp2d_rannacher_turek",2);
    }
#endif

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
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=1;
      const int q=2*k;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      //typedef Dune::PDELab::P1LocalFiniteElementMap<DF,R,GridType::dimension> FEM;
      //FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,GridType::dimension> FEM0;
      FEM0 fem0(Dune::GeometryType::simplex);

      dispatcher(problem,gv,fem,fem0,"UG3d_P1",q);
    }
#endif

#if HAVE_ALUGRID
    // ALU Pk 3D test
    if (false)
    {
      typedef Dune::ALUSimplexGrid<3,3> GridType;

      // read gmsh file
      Dune::GridFactory<GridType> factory(helper.getCommunicator());
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;
      if (helper.rank()==0)
        {
          //Dune::GmshReader<GridType>::read(factory,"grids/cube1045.msh",true,false);
          //Dune::GmshReader<GridType>::read(factory,"grids/rad.msh",true,false);
          Dune::GmshReader<GridType>::read(factory,"grids/plastik-doeddel.msh",true,false);
        }
      GridType *grid=factory.createGrid();
      grid->loadBalance();
      std::cout << " after load balance /" << helper.rank() << "/ " << grid->size(0) << std::endl;
      //grid->globalRefine(1);
      std::cout << " after refinement /" << helper.rank() << "/ " << grid->size(0) << std::endl;

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid->leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=1;
      const int q=2*k;
      typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      //typedef Dune::PDELab::P1LocalFiniteElementMap<DF,R,GridType::dimension> FEM;
      //FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,GridType::dimension> FEM0;
      FEM0 fem0(Dune::GeometryType::simplex);

      dispatcher(problem,gv,fem,fem0,"ALU3d_P1_PlastkDoeddel",q);
    }

    // ALU Q1 3D test
    if (false)
    {
      // make grid
      typedef Dune::ALUCubeGrid<3,3> GridType;
      Dune::GridFactory<GridType> factory(helper.getCommunicator());

      if (helper.rank()==0)
        {
          Dune::FieldVector< double, 3 > pos;
          pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
          pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    factory.insertVertex(pos);
          pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
          pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    factory.insertVertex(pos);
          pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
          pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    factory.insertVertex(pos);
          pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);
          pos[0] = 1;  pos[1] = 1;  pos[2] = 1;    factory.insertVertex(pos);

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
      const GV& gv=grid->leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,R,GridType::dimension> FEM;
      FEM fem;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,GridType::dimension> FEM0;
      FEM0 fem0(Dune::GeometryType::cube);
 
      dispatcher(problem,gv,fem,fem0,"ALU3d_Q1",2);
    }
#endif

#if HAVE_ALBERTA
    if (true)
    {
      typedef AlbertaUnitSquare GridType;
      GridType grid;
      grid.globalRefine(8);

      // get view
      typedef GridType::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GridType::ctype DF;
      typedef double R;
      const int k=4;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,GridType::dimension> FEM0;
      FEM0 fem0(Dune::GeometryType::simplex);

      dispatcher(problem,gv,fem,fem0,"Alberta_Pk",q);
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
