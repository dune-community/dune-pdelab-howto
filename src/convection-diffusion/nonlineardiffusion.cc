// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Nonlinear diffusion equation solved in sequential and parallel (overlapping and nonoverlapping) using linear elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/convectiondiffusion.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include<dune/pdelab/instationary/onestep.hh>

#include"../utility/gridexamples.hh"

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
    return -1.0;
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
  template<typename I>
  bool isDirichlet(
				   const I & intersection,               /*@\label{bcp:name}@*/
				   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
				   ) const
  {
    return true;  // Dirichlet b.c. everywhere
    /*
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );

    if( xg[0]<1E-6 || xg[0]>1-1E-6 )
      return true;  // Dirichlet b.c.
    else
      return false; // Neumann b.c.
    */
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

  //! set time for subsequent evaluation
  void setTime (RF t)
  {}
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

  // <<<2>>> Make grid function space
  const int degree=2;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::BCTypeParam_CD<Param> bctype(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints( bctype, gfs, cg );
  std::cout << "constrained dofs=" << cg.size()
            << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP;
  LOP lop(param,8);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(27); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO;
  GO go(gfs,cg,gfs,cg,lop,mbe);

  // <<<5>>> Compute affine shift
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,x);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write("nonlineardiffusion_initial_guess_Q2",Dune::VTK::ascii);
  }

  // <<<6>>> Make a linear solver
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls(false);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,V> newton(go,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<7b>>> example for solving a linear problem
  //Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> lp(go,x,ls,1E-10);
  //lp.apply();

  // <<<8>>> graphical output
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,x);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write("nonlineardiffusion_sequential_Q2",Dune::VTK::ascii);
  }
}

// a parallel variant for nonoverlapping grids
template<class GV>
void parallel_nonoverlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  int rank = gv.comm().rank();

  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints<GV> CON;
  CON con(gv);
  // <<<2>>> Make grid function space
  const int degree=1;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);

  con.compute_ghosts(gfs); // con stores indices of ghost dofs

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::BCTypeParam_CD<Param> bctype(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints( bctype, gfs, cg );
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size()
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP;
  LOP lop(param,2);
  typedef Dune::PDELab::ISTLMatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,C,C,true> GO;
  GO go(gfs,cg,gfs,cg,lop);

  // <<<5>>> Compute affine shift
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<6>>> Make a linear solver
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS;
  LS ls (go,5000,3,2);
  //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
  //LS ls(gfs,5000,1);

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,V> newton(go,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_nonoverlapping_Q1",Dune::VTK::ascii);
}

// a parallel variant for overlapping grids
template<class GV>
void parallel_overlapping_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  int rank = gv.comm().rank();

  // <<<2>>> Make grid function space
  const int degree=1;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> define problem parameters
  typedef ConvectionDiffusionProblem<GV,Real> Param;
  Param param;
  Dune::PDELab::BCTypeParam_CD<Param> bctype(gv,param);
  typedef Dune::PDELab::DirichletBoundaryCondition_CD<Param> G;
  G g(gv,param);

  // <<<3>>> Compute constrained space
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  Dune::PDELab::constraints( bctype, gfs, cg );
  if (rank==0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cg.size()
                         << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::ConvectionDiffusion<Param> LOP;
  LOP lop(param,2);
  typedef Dune::PDELab::ISTLMatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO;
  GO go(gfs,cg,gfs,cg,lop);

  // <<<4>>> Compute affine shift
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // <<<6>>> Make a linear solver
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> LS;
  LS ls(gfs,cg,5000,1);
#else
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,C> LS;
  LS ls(gfs,cg,5000,10,1);
#endif

  // <<<7>>> solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,V> newton(go,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.apply();

  // <<<8>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("nonlineardiffusion_overlapping_Q1",Dune::VTK::ascii);
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

	if (argc!=2)
	  {
		if(helper.rank()==0)
		  std::cout << "usage: ./nonlineardiffusion <level>" << std::endl;
		return 1;
	  }

	int level;
	sscanf(argv[1],"%d",&level);

    // sequential version
    if (helper.size()==1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(16));
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      sequential(gv);
    }

    // nonoverlapping version
    if (helper.size()>1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(16));
      std::bitset<2> periodic(false);
      int overlap=0; // needs overlap 0 because overlap elements are not assembled anyway
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      parallel_nonoverlapping_Q1(gv);
    }

    // overlapping version
    if (helper.size()>1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(16));
      std::bitset<2> periodic(false);
      int overlap=2;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      parallel_overlapping_Q1(gv);
    }

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
