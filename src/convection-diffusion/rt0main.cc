// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Problems A-F using lowest order Raviart-Thomas elements (sequential)
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<dune/common/parallel/mpihelper.hh>
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
#include<dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/raviartthomas0.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
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

class BCTypeParam_Dummy
  : public Dune::PDELab::DirichletConstraintsParameters /*@\label{bcp:base}@*/
{
public:

  template<typename I>
  bool isDirichlet(
				   const I & intersection,   /*@\label{bcp:name}@*/
				   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
				   ) const
  {
	return false;
  }
};

//===============================================================
// Problem setup and solution
//===============================================================

template<typename BCType, typename GType, typename KType, typename A0Type, typename FType, typename VType,
         typename GV, typename PFEM, typename VFEM>
void driver (BCType& bctype, GType& g, KType& k, A0Type& a0, FType& f, VType& v,
             const GV& gv, const PFEM& pfem, const VFEM& vfem, std::string filename)
{
  // types and constants
  typedef double R;

  // make a grid function space
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,PFEM,Dune::PDELab::NoConstraints,
                                          VBE> P0GFS;
  P0GFS p0gfs(gv,pfem);
  typedef Dune::PDELab::GridFunctionSpace<GV,VFEM,Dune::PDELab::RT0Constraints,
                                          VBE> RT0GFS;
  RT0GFS rt0gfs(gv,vfem);
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VBE,
    Dune::PDELab::LexicographicOrderingTag,
    RT0GFS,P0GFS> MGFS;
  MGFS mgfs(rt0gfs,p0gfs); // the mixed grid function space

  // construct a composite boundary condition type function
  BCTypeParam_Dummy d;
  typedef Dune::PDELab::CompositeConstraintsParameters<BCType,BCTypeParam_Dummy> BCT;
  BCT bct(bctype,d);

  // constraints
  typedef typename MGFS::template ConstraintsContainer<R>::Type T;
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

  // make grid operator
  typedef Dune::PDELab::DiffusionMixed<KType,A0Type,FType,BCType,GType> LOP;
  LOP lop(k,a0,f,bctype,g,4,2);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeros per row can be cross-checked by patternStatistics().
  typedef Dune::PDELab::GridOperator<MGFS,MGFS,LOP,MBE,R,R,R,T,T> GO;
  GO go(mgfs,t,mgfs,t,lop,mbe);

  // represent operator as a matrix
  typedef typename GO::Jacobian M;
  M m(go);
  std::cout << m.patternStatistics() << std::endl;
  m = 0.0;
  go.jacobian(x,m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // set up solver
  typedef typename M::BaseT ISTLM;
#if HAVE_SUPERLU
  Dune::SuperLU<ISTLM> solver(Dune::PDELab::istl::raw(m), true);
  Dune::InverseOperatorResult stat;

  X r(mgfs,0.0);
  go.residual(x,r);
  X z(mgfs,0.0);
  solver.apply(z,r,stat);
  x -= z;
#else
#error No superLU support, please install and configure it.
#endif

  // select subspaces
  typedef Dune::PDELab::GridFunctionSubSpace
    <MGFS,Dune::TypeTree::TreePath<0> > VSUB;

  VSUB vsub(mgfs);                   // velocity subspace

  typedef Dune::PDELab::GridFunctionSubSpace
      <MGFS,Dune::TypeTree::TreePath<1> > PSUB;
  PSUB psub(mgfs);                   // pressure subspace

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionPiola<VSUB,X> RT0DGF;
  RT0DGF rt0dgf(vsub,x);
  typedef Dune::PDELab::DiscreteGridFunction<PSUB,X> P0DGF;
  P0DGF p0dgf(psub,x);

  // output grid function with VTKWriter
  //Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1); // plot result
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P0DGF>(p0dgf,"pressure"));
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>(rt0dgf,"velocity"));
  vtkwriter.write(filename,Dune::VTK::ascii);
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
      BCTypeParam_A bctype;
      G_A<GV,RF> g(gv);
      K_A<GV,RF> k(gv);
      A0_A<GV,RF> a0(gv);
      F_A<GV,RF> f(gv);
      V_A<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==B)
    {
      BCTypeParam_B bctype;
      G_B<GV,RF> g(gv);
      K_B<GV,RF> k(gv);
      A0_B<GV,RF> a0(gv);
      F_B<GV,RF> f(gv);
      V_B<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==C)
    {
      BCTypeParam_C bctype;
      G_C<GV,RF> g(gv);
      K_C<GV,RF> k(gv);
      A0_C<GV,RF> a0(gv);
      F_C<GV,RF> f(gv);
      V_C<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==D)
    {
      Dune::FieldVector<RF,GV::Grid::dimension> correlation_length;
      correlation_length = 1.0/64.0;

      BCTypeParam_D bctype;
      G_D<GV,RF> g(gv);
      K_D<GV,RF> k(gv,correlation_length,1.0,0.0,5000,-1083);
      A0_D<GV,RF> a0(gv);
      F_D<GV,RF> f(gv);
      V_D<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==E)
    {
      BCTypeParam_E bctype;
      G_E<GV,RF> g(gv);
      K_E<GV,RF> k(gv);
      A0_E<GV,RF> a0(gv);
      F_E<GV,RF> f(gv);
      V_E<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
    }
  if (problem==F)
    {
      BCTypeParam_F bctype;
      G_F<GV,RF> g(gv);
      K_F<GV,RF> k(gv);
      A0_F<GV,RF> a0(gv);
      F_F<GV,RF> f(gv);
      V_F<GV,RF> v(gv);
      driver(bctype,g,k,a0,f,v,gv,pfem,vfem,filename);
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
      const int dim = 2;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(7);

      // instantiate finite element maps
      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::YaspGrid<dim>::LeafGridView GV;

      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,R,0,Dune::GeometryType::cube> RT0FEM;
      RT0FEM rt0fem(grid.leafGridView());

      dispatcher(problem,grid.leafGridView(),p0fem,rt0fem,"Yasp2d_rt0q");
    }

    // YaspGrid 3D test
    if (true)
    {
      const int dim = 3;
      Dune::FieldVector<double,dim> L(1.0);
      Dune::array<int,dim> N(Dune::fill_array<int,dim>(1));
      std::bitset<dim> B(false);
      Dune::YaspGrid<dim> grid(L,N,B,0);
      grid.globalRefine(4);

      // instantiate finite element maps
      typedef Dune::YaspGrid<dim>::ctype DF;
      typedef Dune::YaspGrid<dim>::LeafGridView GV;

      typedef double R;
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,R,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube,dim));

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,R,0,Dune::GeometryType::cube> RT0FEM;

      RT0FEM rt0fem(grid.leafGridView());

      dispatcher(problem,grid.leafGridView(),p0fem,rt0fem,"Yasp3d_rt0q");
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

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,R,0,Dune::GeometryType::simplex> RT0FEM;

      RT0FEM rt0fem(grid.leafGridView());

      dispatcher(problem,grid.leafGridView(),p0fem,rt0fem,"ALU2d_rt0");
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

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,R,0,Dune::GeometryType::simplex> RT0FEM;

      RT0FEM rt0fem(grid.leafGridView());

      dispatcher(problem,grid.leafGridView(),p0fem,rt0fem,"UG2d_rt0");
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

      typedef Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,R,0,Dune::GeometryType::simplex> RT0FEM;

      RT0FEM rt0fem(grid.leafGridView());

      dispatcher(problem,grid.leafGridView(),p0fem,rt0fem,"Alberta2d_rt0");
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
