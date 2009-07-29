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
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
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

#include"problemD.hh"

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template<typename BType, typename GType, typename KType, typename A0Type, typename FType, typename JType,
         typename GV, typename FEM> 
void driver (const BType& b, const GType& g, 
             const KType& k, const A0Type& a0, const FType& f, const JType& j,
             const GV& gv, const FEM& fem, std::string filename, int intorder=1)
{
  // constants and types and global variables
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::ConformingDirichletConstraints,
    Dune::PDELab::ISTLVectorBackend<1>,
    Dune::PDELab::SimpleGridFunctionStaticSize
    > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<R>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // make grid function operator
  typedef Dune::PDELab::Diffusion<KType,A0Type,FType,BType,JType> LOP; 
  LOP lop(k,a0,f,b,j,intorder);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
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
  POP pop(gfs,m,phelper);
  typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
  PSP psp(gfs,phelper);
  typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,V> PRICH;
  PRICH prich(gfs,phelper);
  int verbose;
  if (gv.comm().rank()==0) verbose=1; else verbose=0;
  Dune::CGSolver<V> solver(pop,psp,prich,1E-8,40000,verbose);
  Dune::InverseOperatorResult stat;  
  solver.apply(z,r,stat);
  x -= z;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  //vtkwriter.pwrite(filename.c_str(),"vtk","",Dune::VTKOptions::binaryappended);
  vtkwriter.write(filename,Dune::VTKOptions::ascii);
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
    
#if HAVE_MPI
    // Q1, 2d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(4);
      Dune::FieldVector<bool,2> B(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
      grid.globalRefine(4);
      
      typedef Dune::YaspGrid<2>::ctype DF;
      typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;

      Dune::FieldVector<double,2> correlation_length;
      correlation_length = 1.0/64.0;
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();

      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_yasp2d_Q1");
    }
#endif

#if HAVE_MPI
    // Q1, 3d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(1);
      Dune::FieldVector<bool,3> B(false);
      int overlap=0;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
      grid.globalRefine(2);

      typedef Dune::YaspGrid<3>::ctype DF;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;

      Dune::FieldVector<double,3> correlation_length;
      correlation_length = 1.0/64.0;
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafView();

      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_yasp3d_Q1");
    }
#endif

#if HAVE_ALUGRID
    // ALU Pk 3D test
    if (true)
    {
      typedef Dune::ALUSimplexGrid<3,3> GridType;

      // read gmsh file
      Dune::GridFactory<GridType> factory(helper.getCommunicator());
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;
      if (helper.rank()==0)
        Dune::GmshReader<GridType>::read(factory,"grids/cube1045.msh",true,false);
      GridType *grid=factory.createGrid();
      grid->loadBalance();
      std::cout << " after load balance /" << helper.rank() << "/ " << grid->size(0) << std::endl;
      grid->globalRefine(1);
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

      Dune::FieldVector<double,3> correlation_length;
      correlation_length = 1.0/64.0;

      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_ALU3d_P1",q);
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

      Dune::FieldVector<double,3> correlation_length;
      correlation_length = 1.0/64.0;

      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_ALU3d_Q1");
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
