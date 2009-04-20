// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/rt02dfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>

#include"gridexamples.hh"
#include"rt0constraints.hh"
#include"poissonhdivconforming.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f         in \Omega, 
//                    u = g         on \partial\Omega_D
//  -\nabla u \cdot \nu = v\cdot\nu on \partial\Omega_N
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
  }
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<GV,int,1,
                                                                             Dune::FieldVector<int,1> >,
                                                  B<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B<GV> > BaseT;

  B (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      xg = ig.geometry().global(x);

    if (xg[1]<1E-6 || xg[1]>1.0-1E-6)
      {
        y = 0; // Neumann
        return;
      }
    if (xg[0]>1.0-1E-6 && xg[1]>0.5+1E-6)
      {
        y = 0; // Neumann
        return;
      }
    y = 1; // Dirichlet
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

template<typename GV>
class Dummy
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<GV,int,1,
                                                                             Dune::FieldVector<int,1> >,
                                                  Dummy<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,Dummy<GV> > BaseT;

  Dummy (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {}

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
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

template<typename GV, typename RF>
class V 
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2>,
													  V<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,2> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V<GV,RF> > BaseT;

  V (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {  
    if (x[1]<1E-6 || x[1]>1.0-1E-6)
      {
        y[0] = 0; y[1] = 0;
        return;
      }
    if (x[0]>1.0-1E-6 && x[1]>0.5+1E-6)
      {
        y[0] = -5.0; y[1] = 0.0;
        return;
      }
    y = 0.0;
  }
};

//===============================================================
// Problem setup and solution 
//===============================================================

template<class GV> 
void testrt0 (const GV& gv, std::string filename)
{
  // types and constants
  typedef typename GV::Grid::ctype DF;
  const int dim = GV::dimension;
  typedef double R;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
  P0FEM p0fem(Dune::GeometryType::simplex);
  typedef Dune::PDELab::RT02DLocalFiniteElementMap<GV,DF,R> RT0FEM;
  RT0FEM rt0fem(gv);
  
  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,
    Dune::PDELab::ISTLVectorBackend<1> > P0GFS; 
  P0GFS p0gfs(gv,p0fem);
  typedef Dune::PDELab::GridFunctionSpace<GV,RT0FEM,RT0Constraints,
    Dune::PDELab::ISTLVectorBackend<1> > RT0GFS; 
  RT0GFS rt0gfs(gv,rt0fem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::GridFunctionSpaceLexicographicMapper,
	  RT0GFS,P0GFS> MGFS;              
  MGFS mgfs(rt0gfs,p0gfs); // the mixed grid function space

  // construct a composite boundary condition type function
  typedef B<GV> BType;
  BType b(gv);
  typedef Dummy<GV> DType;
  DType d(gv);
  typedef Dune::PDELab::CompositeGridFunction<BType,DType> BCT;
  BCT bct(b,d);

  // constraints 
  typedef typename RT0GFS::template ConstraintsContainer<R>::Type T;
  T t;                               // container for transformation    
  Dune::PDELab::constraints(bct,mgfs,t); // fill container

  // construct a composite grid function
  typedef V<GV,R> VType;
  VType v(gv);
  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rv(v);
  typedef G<GV,R> GType;
  GType g(gv);
  typedef Dune::PDELab::CompositeGridFunction<RVType,GType> UType;
  UType u(rv,g);

  // make coefficent Vectors
  typedef typename MGFS::template VectorContainer<R>::Type X;
  X x(mgfs,0.0);

  // do interpolation
  Dune::PDELab::interpolate(u,mgfs,x);
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x);  // clear interior

  // make grid operator space
  typedef F<GV,R> FType;
  FType f(gv);
  typedef Dune::PDELab::PoissonHDivConforming<FType,BType,GType> LOP; 
  LOP lop(f,b,g);
  typedef Dune::PDELab::GridOperatorSpace<MGFS,MGFS,
    LOP,T,T,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(mgfs,t,mgfs,t,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  M m(gos);
  m = 0.0;
  gos.jacobian(x,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // set up solver
  typedef typename M::BaseT ISTLM;
  Dune::SuperLU<ISTLM> solver(m, true);
  Dune::InverseOperatorResult stat;

  X r(mgfs,0.0);
  gos.residual(x,r);
  X z(mgfs,0.0);
  solver.apply(z,r,stat);
  x -= z;

  // select subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS,0> VSUB;
  VSUB vsub(mgfs);                   // velocity subspace
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS,1> PSUB;
  PSUB psub(mgfs);                   // pressure subspace

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionPiola<VSUB,X> RT0DGF;
  RT0DGF rt0dgf(vsub,x);
  typedef Dune::PDELab::DiscreteGridFunction<PSUB,X> P0DGF;
  P0DGF p0dgf(psub,x);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,4); // plot result
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P0DGF>(p0dgf,"pressure"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>(rt0dgf,"velocity"));
  vtkwriter.write(filename,Dune::VTKOptions::ascii);
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

#if HAVE_ALUGRID
    {
      ALUUnitSquare grid;
      grid.globalRefine(6);
      testrt0(grid.leafView(),"rt0_ALU");
    }
#endif

#if HAVE_UG
    {
      UGUnitSquare grid;
      grid.globalRefine(6);
      testrt0(grid.leafView(),"rt0_UG");
    }
#endif

#if HAVE_ALBERTA
    {
      AlbertaUnitSquare grid;
      grid.globalRefine(12);
      testrt0(grid.leafView(),"rt0_Alberta");
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
