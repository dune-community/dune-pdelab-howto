// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Solve Problems A-F using lowest order Raviart-Thomas elements (sequential)
*/
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
#include<dune/pdelab/finiteelementmap/rt0qfem.hh>
#include<dune/pdelab/finiteelementmap/rt0constraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/localoperator/diffusionmixed.hh>

#include"../utility/gridexamples.hh"
#include"problemA.hh"
#include"problemB.hh"
#include"problemC.hh"
#include"problemD.hh"
#include"problemE.hh"
#include"problemF.hh"

//===============================================================
// dummy boundary condition function for the pressure component
//===============================================================

template<typename GV>
class Dummy
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
               BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,
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

//===============================================================
// Problem setup and solution 
//===============================================================

template<typename BType, typename GType, typename KType, typename A0Type, typename FType, typename VType, 
         typename GV, typename PFEM, typename VFEM> 
void driver (BType& b, GType& g, KType& k, A0Type& a0, FType& f, VType& v,
             const GV& gv, const PFEM& pfem, const VFEM& vfem, std::string filename)
{
  // types and constants
  typedef typename GV::Grid::ctype DF;
  typedef double R;

  // make a grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,PFEM,Dune::PDELab::NoConstraints,
    Dune::PDELab::ISTLVectorBackend<1> > P0GFS; 
  P0GFS p0gfs(gv,pfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,VFEM,Dune::PDELab::RT0Constraints,
    Dune::PDELab::ISTLVectorBackend<1> > RT0GFS; 
  RT0GFS rt0gfs(gv,vfem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::GridFunctionSpaceLexicographicMapper,
	  RT0GFS,P0GFS> MGFS;              
  MGFS mgfs(rt0gfs,p0gfs); // the mixed grid function space

  // construct a composite boundary condition type function
  typedef Dummy<GV> DType;
  DType d(gv);
  typedef Dune::PDELab::CompositeGridFunction<BType,DType> BCT;
  BCT bct(b,d);

  // constraints 
  typedef typename RT0GFS::template ConstraintsContainer<R>::Type T;
  T t;                               // container for transformation    
  Dune::PDELab::constraints(bct,mgfs,t); // fill container

  // construct a composite grid function
  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rv(v);
  typedef Dune::PDELab::CompositeGridFunction<RVType,GType> UType;
  UType u(rv,g);

  // make coefficent Vectors
  typedef typename Dune::PDELab::BackendVectorSelector<MGFS,R>::Type X;
  X x(mgfs,0.0);

  // do interpolation
  Dune::PDELab::interpolate(u,mgfs,x);
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x);  // clear interior

  // make grid operator space
  typedef Dune::PDELab::DiffusionMixed<KType,A0Type,FType,BType,GType> LOP; 
  LOP lop(k,a0,f,b,g,4,2);
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
#if HAVE_SUPERLU
  Dune::SuperLU<ISTLM> solver(m, true);
  Dune::InverseOperatorResult stat;

  X r(mgfs,0.0);
  gos.residual(x,r);
  X z(mgfs,0.0);
  solver.apply(z,r,stat);
  x -= z;
#else
#error No superLU support, please install and configure it.
#endif

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
  //Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1); // plot result
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<P0DGF>(p0dgf,"pressure"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>(rt0dgf,"velocity"));
  vtkwriter.write(filename,Dune::VTKOptions::ascii);
}

template<typename GV, typename PFEM, typename VFEM> 
void dispatcher (std::string problem, const GV& gv, const PFEM& pfem, const VFEM& vfem, std::string gridname)
{
  std::string A("A"), B("B"), C("C"), D("D"), E("E"), F("F");
  std::string filename(""), underscore("_");
  filename = problem+underscore+gridname;

  typedef double RF;

  if (problem==A) 
    {
      B_A<GV> b(gv); 
      G_A<GV,RF> g(gv);
      K_A<GV,RF> k(gv);
      A0_A<GV,RF> a0(gv);
      F_A<GV,RF> f(gv);
      V_A<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==B) 
    {
      B_B<GV> b(gv); 
      G_B<GV,RF> g(gv);
      K_B<GV,RF> k(gv);
      A0_B<GV,RF> a0(gv);
      F_B<GV,RF> f(gv);
      V_B<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==C) 
    {
      B_C<GV> b(gv); 
      G_C<GV,RF> g(gv);
      K_C<GV,RF> k(gv);
      A0_C<GV,RF> a0(gv);
      F_C<GV,RF> f(gv);
      V_C<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==D) 
    {
      Dune::FieldVector<RF,GV::Grid::dimension> correlation_length;
      correlation_length = 1.0/64.0;

      B_D<GV> b(gv); 
      G_D<GV,RF> g(gv);
      K_D<GV,RF> k(gv,correlation_length,1.0,0.0,5000,-1083);
      A0_D<GV,RF> a0(gv);
      F_D<GV,RF> f(gv);
      V_D<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==E) 
    {
      B_E<GV> b(gv); 
      G_E<GV,RF> g(gv);
      K_E<GV,RF> k(gv);
      A0_E<GV,RF> a0(gv);
      F_E<GV,RF> f(gv);
      V_E<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==F) 
    {
      B_F<GV> b(gv); 
      G_F<GV,RF> g(gv);
      K_F<GV,RF> k(gv);
      A0_F<GV,RF> a0(gv);
      F_F<GV,RF> f(gv);
      V_F<GV,RF> v(gv);
      driver(b,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    std::string problem="A";

    // YaspGrid 2D test
    if (true)
    {  
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(7);

      // instantiate finite element maps
      typedef Dune::YaspGrid<2>::ctype DF;
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const int dim = 2;
      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
      typedef Dune::PDELab::RT0QLocalFiniteElementMap<GV,DF,R,dim> RT0FEM;
      RT0FEM rt0fem(grid.leafView());

      dispatcher(problem,grid.leafView(),p0fem,rt0fem,"Yasp2d_rt0q");
    }

    // YaspGrid 3D test
    if (true)
    {
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(1);
      Dune::FieldVector<bool,3> B(false);
      Dune::YaspGrid<3> grid(L,N,B,0);
      grid.globalRefine(4);
 
      // instantiate finite element maps
      typedef Dune::YaspGrid<3>::ctype DF;
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const int dim = 3;
      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
      typedef Dune::PDELab::RT0QLocalFiniteElementMap<GV,DF,R,dim> RT0FEM;
      RT0FEM rt0fem(grid.leafView());

      dispatcher(problem,grid.leafView(),p0fem,rt0fem,"Yasp3d_rt0q");
    }

#if HAVE_ALUGRID
    if (true)
    {
      ALUUnitSquare grid;
      grid.globalRefine(7);

      // instantiate finite element maps
      typedef ALUUnitSquare::ctype DF;
      typedef ALUUnitSquare::LeafGridView GV;
      const int dim = 2;
      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
      typedef Dune::PDELab::RT02DLocalFiniteElementMap<GV,DF,R> RT0FEM;
      RT0FEM rt0fem(grid.leafView());

      dispatcher(problem,grid.leafView(),p0fem,rt0fem,"ALU2d_rt0");
    }
#endif

#if HAVE_UG
    if (true)
    {
      UGUnitSquare grid;
      grid.globalRefine(6);

      // instantiate finite element maps
      typedef UGUnitSquare::ctype DF;
      typedef UGUnitSquare::LeafGridView GV;
      const int dim = 2;
      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
      typedef Dune::PDELab::RT02DLocalFiniteElementMap<GV,DF,R> RT0FEM;
      RT0FEM rt0fem(grid.leafView());

      dispatcher(problem,grid.leafView(),p0fem,rt0fem,"UG2d_rt0");
    }
#endif

#if HAVE_ALBERTA
    if (true)
    {
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // instantiate finite element maps
      typedef AlbertaUnitSquare::ctype DF;
      typedef AlbertaUnitSquare::LeafGridView GV;
      const int dim = 2;
      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
      typedef Dune::PDELab::RT02DLocalFiniteElementMap<GV,DF,R> RT0FEM;
      RT0FEM rt0fem(grid.leafView());

      dispatcher(problem,grid.leafView(),p0fem,rt0fem,"Alberta2d_rt0");
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
