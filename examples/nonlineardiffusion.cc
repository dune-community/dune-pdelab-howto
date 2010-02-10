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
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

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
#include<dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/diffusion.hh>
#include<dune/pdelab/localoperator/convectiondiffusion.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include<dune/pdelab/instationary/onestep.hh>

#include"gridexamples.hh"

//==============================================================================
// Parameter class for the convection diffusion problem
//==============================================================================

//! base class for parameter class
template<typename GV, typename RF>
class ConvectionDiffusionProblem : 
  public Dune::PDELab::ConvectionDiffusionParameterInterface<Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>, 
                                                             ConvectionDiffusionProblem<GV,RF> >
{
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! source/reaction term
  typename Traits::RangeFieldType 
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return 0.0;
  }

  //! nonlinearity under gradient
  typename Traits::RangeFieldType 
  w (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return u;
  }

  //! nonlinear scaling of diffusion tensor
  typename Traits::RangeFieldType 
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    return 1.0;
  }

  //! tensor permeability
  typename Traits::PermTensorType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType kabs;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        kabs[i][j] = (i==j) ? 1 : 0;
    return kabs;
  }

  //! nonlinear flux vector
  typename Traits::RangeType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
     typename Traits::RangeFieldType u) const
  {
    typename Traits::RangeType flux;
    flux[0] = 50 * 1.00 * u*u;
    flux[1] = 50 * 0.25 * u*u;
    return flux;
  }

  //! boundary condition type function
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  int
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::RangeType global = is.geometry().global(x);
    return 1; 
    if (global[0]<1E-6 || global[0]>1-1E-6)
      return 1;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType global = e.geometry().global(x);
    if (global[0]>1E-12 && global[0]<1.0-1E-12 &&           // this is needed for parallel
        global[1]>1E-12 && global[1]<1.0-1E-12) return 0.0; // overlapping version !
    if (global[0]<1E-6 && global[1]>0.25 && global[1]<0.75)
      return 1.0;
    else
      return 0.0;
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType 
  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType u) const
  {
    return 0.0;
  }
};


//===============================================================
// Some variants to solve the nonlinear diffusion problem
//===============================================================

// a sequential variant
template<class GV>
void sequential (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_CD<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints(b,gfs,cg);
  std::cout << "constrained dofs=" << cg.size() 
            << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP; 
  LOP lop(param,6);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
#ifdef HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls(false);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<7b>>> example for solving a linear problem
  //Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,V> lp(gos,x,ls,1E-10);
  //lp.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_sequential_Q1",Dune::VTKOptions::ascii);
}

// a parallel variant for nonoverlapping grids
template<class GV>
void parallel_nonoverlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  int rank = gv.comm().rank();

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem,con);
  con.compute_ghosts(gfs); // con stores indices of ghost dofs

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_CD<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints(b,gfs,cg);
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size() 
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP; 
  LOP lop(param,2);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE,true> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
  LS ls(gfs,5000,1);

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_nonoverlapping_Q1",Dune::VTKOptions::ascii);
}

// a parallel variant for overlapping grids
template<class GV>
void parallel_overlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  int rank = gv.comm().rank();

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  typedef Dune::PDELab::BoundaryConditionType_CD<Param> B;
  B b(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);
  
  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints(b,gfs,cg);
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size() 
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Compute affine shift
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP; 
  LOP lop(param,2);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,C,C,MBE> GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // <<<6>>> Make a linear solver 
#ifdef HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> LS;
  LS ls(gfs,cg,5000,2);
#else
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,C> LS;
  LS ls(gfs,cg,5000,10,2);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,V> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_overlapping_Q1",Dune::VTKOptions::ascii);
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

    // sequential version
    if (1 && helper.size()==1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(8);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(3);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      sequential(gv);
    }

#if HAVE_MPI
    // nonoverlapping version
    if (1 && helper.size()>1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(16);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=0; // needs overlap 0 because overlap elements are not assembled anyway
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      grid.globalRefine(2);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      parallel_nonoverlapping_Q1(gv);
    }

    // overlapping version
    if (1 && helper.size()>1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(8);
      Dune::FieldVector<bool,2> periodic(false);
      int overlap=2;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      grid.globalRefine(3);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      parallel_overlapping_Q1(gv);
    }
#endif

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
