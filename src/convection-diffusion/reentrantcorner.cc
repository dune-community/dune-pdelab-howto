// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief Solve Poisson problem in reentrant corner domain (sequential)
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
#include<dune/grid/common/gridfactory.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#include<dune/grid/io/file/dgfparser/dgfug.hh>
#include<dune/grid/uggrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/dgfparser/dgfyasp.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/poisson.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"../utility/gridexamples.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
// In this case we actually want to solve the laplace equation,
// i.e.
//                 f = 0
//  \partial\Omega_N = \emptyset
//                 j = arbitrary
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
    y = 0.0;
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
    //typename Traits::DomainFieldType rho = std::sqrt(x[0]*x[0]+x[1]*x[1]);
    typename Traits::DomainFieldType theta = std::atan2(x[1], x[0]);
    if(theta < 0.0) theta += 2*M_PI;
    // how much pacman opens its mouth
    const typename Traits::DomainFieldType opening = 10.0/180.0*M_PI;
    if(theta > 2*M_PI - opening) {
      y = 0.0;
    }
    else {
      y = std::sin(theta/(2*M_PI-opening)*M_PI);
    }
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
    y = 0.0;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON>
void poisson_driver (const GV& gv, const FEM& fem, std::string filename, int k)
{
  const int q=2*k; // integration order

  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    VBE > GFS;
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();

  typedef Dune::PDELab::DirichletConstraintsParameters BCTypeParam;
  BCTypeParam bctype;
  Dune::PDELab::constraints(bctype,gfs,cg);

  // make grid function operator
  typedef F<GV,R> FType;
  FType f(gv);
  typedef J<GV,R> JType;
  JType j(gv);
  typedef Dune::PDELab::Poisson<FType,BCTypeParam,JType> LOP;
  LOP lop(f,bctype,j,q);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

  // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
  // Unstructured grid with triangular elements where one node can be part of up to 9 cells.
  MBE mbe(60);

  typedef Dune::PDELab::GridOperator<GFS,GFS,
                                     LOP,MBE,R,R,R,C,C> GO;
  GO go(gfs,cg,gfs,cg,lop,mbe);

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  typedef G<GV,R> GType;
  GType g(gv);
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  /* create solver yourself (for example to look at global stiffness matrix)
  // represent operator as a matrix
  typedef typename GO::template MatrixContainer<R>::Type M;
  M m(go);
  m = 0.0;
  go.jacobian(x0,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // evaluate residual w.r.t initial guess
  V r(gfs);
  r = 0.0;
  go.residual(x0,r);

  // make ISTL solver
  Dune::MatrixAdapter<M,V,V> opa(m);
  typedef Dune::PDELab::OnTheFlyOperator<V,V,GO> ISTLOnTheFlyOperator;
  ISTLOnTheFlyOperator opb(go);
  Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
  Dune::SeqILU0<M,V,V> ilu0(m,1.0);
  Dune::Richardson<V,V> richardson(1.0);

  //   typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<M,
  //     Dune::Amg::FirstDiagonal> > Criterion;
  //   typedef Dune::SeqSSOR<M,V,V> Smoother;
  //   typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
  //   SmootherArgs smootherArgs;
  //   smootherArgs.iterations = 2;
  //   int maxlevel = 20, coarsenTarget = 100;
  //   Criterion criterion(maxlevel, coarsenTarget);
  //   criterion.setMaxDistance(2);
  //   typedef Dune::Amg::AMG<Dune::MatrixAdapter<M,V,V>,V,Smoother> AMG;
  //   AMG amg(opa,criterion,smootherArgs,1,1);

  Dune::CGSolver<V> solvera(opa,ilu0,1E-10,5000,2);
  Dune::CGSolver<V> solverb(opb,richardson,1E-10,5000,2);
  Dune::InverseOperatorResult stat;

  // solve the jacobian system
  r *= -1.0; // need -residual
  V x(gfs,0.0);
  solvera.apply(x,r,stat);
  x += x0;
  */

  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls (5000,2);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,x0,ls,1e-12);
  slp.apply();

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x0);

  // output grid function with VTKWriter
  {
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,k-1);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(filename,Dune::VTK::ascii);
  }
  {
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    vtkwriter.write(filename+"_grid",Dune::VTK::ascii);
  }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_ALBERTA
    {
      typedef AlbertaReentrantCorner::Grid Grid;
      // make grid
      AlbertaReentrantCorner gridp;
      Grid &grid = *gridp;
      grid.globalRefine(0);

      // get view
      typedef Grid::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,R,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson_driver<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"reentrantcorner_Alberta_Pk_2d",k);
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
