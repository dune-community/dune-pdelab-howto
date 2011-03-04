// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Solve Problems A-F in parallel using cell-centered finite volumes (works on nonoverlapping grids in overlapping mode)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
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

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p0constraints.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/localoperator/diffusionccfv.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#include"../utility/gridexamples.hh"
#include"problemA.hh"
#include"problemB.hh"
#include"problemC.hh"
#include"problemD.hh"
#include"problemE.hh"
#include"problemF.hh"

template<class GV> 
void test (const GV& gv)
{
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim)); // works only for cubes
  
  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::P0ParallelConstraints,
    Dune::PDELab::ISTLVectorBackend<1>,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS; 
  watch.reset();
  GFS gfs(gv,fem);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // make coefficent Vector and initialize it from a function
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
  V x0(gfs);
  x0 = 0.0;
  typedef G_D<GV,RF> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);

  // local operator
  watch.reset();
  typedef k_D<GV,RF> KType;
  Dune::FieldVector<double,dim> correlation_length;
  correlation_length = 1.0/64.0;
  KType k(gv,correlation_length,0.5,0.0,5000,-1083);
  typedef A0_D<GV,RF> A0Type;
  A0Type a0(gv);
  typedef F_D<GV,RF> FType;
  FType f(gv);
  typedef B_D<GV> BType;
  BType b(gv);
  typedef J_D<GV,RF> JType;
  JType j(gv);
  typedef Dune::PDELab::DiffusionCCFV<KType,A0Type,FType,BType,JType,GType> LOP;
  LOP lop(k,a0,f,b,j,g);
  std::cout << "=== local operator setup " <<  watch.elapsed() << " s" << std::endl;

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::constraints(g,gfs,cc,false);

  // grid operator space
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,
    CC,CC,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<RF>::Type M;
  watch.reset();
  M m(gos);
  std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
  m = 0.0;
  watch.reset();
  gos.jacobian(x0,m);
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
  //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  watch.reset();
  gos.residual(x0,r);
  std::cout << "=== residual evaluation " <<  watch.elapsed() << " s" << std::endl;

  // set up parallel solver
  typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;
  PHELPER phelper(gfs);
  typedef Dune::PDELab::OverlappingOperator<CC,M,V,V> POP;
  POP pop(cc,m);
  typedef Dune::PDELab::OverlappingScalarProduct<GFS,V> PSP;
  PSP psp(gfs,phelper);
//   typedef Dune::PDELab::SuperLUSubdomainSolver<GFS,M,V,V> PSUBSOLVE;
//   PSUBSOLVE psubsolve(gfs,m);
  typedef Dune::SeqSSOR<M,V,V> SeqPrec;
  SeqPrec seqprec(m,10,1.0);
  typedef Dune::PDELab::OverlappingWrappedPreconditioner<CC,GFS,SeqPrec> WPREC;
  WPREC  wprec(gfs,seqprec,cc,phelper);
  int verbose;
  if (gv.comm().rank()==0) verbose=1; else verbose=0;
  Dune::CGSolver<V> solver(pop,psp,wprec,1E-8,40000,verbose);
  Dune::InverseOperatorResult stat;  


  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  solver.apply(x,r,stat);
  x += x0;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.pwrite("single_phase_yasp2d_CCFV","vtk","",Dune::VTKOptions::binaryappended);
  //  vtkwriter.write("single_phase_yasp2d_CCFV",Dune::VTKOptions::ascii);
}

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

    // 2D
    if (true)
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(128);
      Dune::FieldVector<bool,2> B(false);
      int overlap=4;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
      //      grid.globalRefine(6);
      
      // solve problem :)
      test(grid.leafView());
    }

    // Q1, 3d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(64);
      Dune::FieldVector<bool,3> B(false);
      int overlap=4;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);

      // solve problem :)
      test(grid.leafView());
    }

    // UG Q1 2D test
#if HAVE_UG
    if (false)
    {
      // make grid 
      UGUnitSquareQ grid(1000);
      grid.globalRefine(6);

      test(grid.leafView());
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
