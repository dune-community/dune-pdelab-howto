// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Laplace equation with cell-centered finite volume method
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

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
// eigen
#if HAVE_EIGEN
#include<dune/pdelab/backend/eigenvectorbackend.hh>
#include<dune/pdelab/backend/eigenmatrixbackend.hh>
#include<dune/pdelab/backend/eigensolverbackend.hh>
#endif
// istl
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"../utility/gridexamples.hh"


// define some grid functions to interpolate from
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

// selecting the boundary condition type
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

    //Dune::FieldVector<typename I::ctype, I::dimension>
    //  xg = intersection.geometry().global( coord );
    return true;  // Dirichlet b.c. on all boundaries
  }
};


template<class GV>
void test (const GV& gv, std::string filename )
{
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim)); // works only for cubes

  // make function space
#ifdef USE_EIGEN
  typedef Dune::PDELab::EigenVectorBackend VBE;
#else
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
#endif
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::NoConstraints,VBE> GFS;
  watch.reset();
  GFS gfs(gv,fem);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  typedef G<GV,RF> GType;
  GType g(gv);

#ifdef USE_EIGEN
  typedef Dune::PDELab::SparseEigenMatrixBackend MBE;
#else
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
#endif

  // make grid function operator
  typedef Dune::PDELab::LaplaceDirichletCCFV<GType> LOP;
  LOP la(g);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     MBE,
                                     RF,RF,RF,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation > GO;
#ifdef USE_EIGEN
  GO go(gfs,gfs,la);
#else
  GO go(gfs,gfs,la,mbe);
#endif

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x0(gfs);
  x0 = 0.0;
  Dune::PDELab::interpolate(g,gfs,x0);


  V x(gfs,0.0);
#ifdef USE_EIGEN
  typedef Dune::PDELab::EigenBackend_BiCGSTAB_Diagonal LS;
  LS ls (5000);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GO> LS;
  LS ls (5000,2);
#endif
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,ls,x,1e-12);
  slp.apply();

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::nonconforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"u"));
  vtkwriter.write( filename.c_str() ,Dune::VTK::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    std::string basename = "laplacedirichletccfv";
#ifdef USE_EIGEN
    basename += "_eigen";
#endif


    // 2D
    if(helper.size()==1){
      // make grid
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(6);

      // solve problem :)
      test(grid.leafGridView(),basename+"_yasp2d");
    }
    else{
      if(helper.rank()==0)
      std::cout<< "This is a sequential program, you cannot run it in parallel." << std::endl;
    }

    // UG Q1 2D test
#if HAVE_UG
    if(helper.size()==1){
      // make grid
      UGUnitSquareQ grid(1000);
      grid.globalRefine(6);

      test(grid.leafGridView(),basename+"_ug2d");
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
