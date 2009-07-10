// -*- tab-width: 4; indent-tabs-mode: nil -*-
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
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#include"gridexamples.hh"
#include"twophaseop.hh"

//==============================================================================
// Problem definition
//==============================================================================

// initial conditions for phase pressures
template<typename GV, typename RF>
class P_l
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_l<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_l<GV,RF> > BaseT;

  P_l (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
 	y = x[0];
  }
};
template<typename GV, typename RF>
class P_g
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_g<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_g<GV,RF> > BaseT;

  P_g (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
 	y = x[1];
  }
};

// parameter class for local operator
template<typename GV, typename RF>
class TwoPhaseParameter 
  : public Dune::PDELab::TwoPhaseParameterInterface<Dune::PDELab::TwoPhaseParameterTraits<GV,RF>, 
                                                    TwoPhaseParameter<GV,RF> >
{
public:
  typedef Dune::PDELab::TwoPhaseParameterTraits<GV,RF> Traits;

  //! constructor
  TwoPhaseParameter ()
  {
    gvector = 0.0; // no gravity
  }

  //! porosity
  typename Traits::RangeFieldType 
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.2;
  }

  //! inverse capillary pressure function
  typename Traits::RangeFieldType 
  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType pc) const
  {
    return 1.0-0.001*pc;
  }
	  
  //! liquid phase relative permeability
  typename Traits::RangeFieldType 
  kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType s_l) const
  {
    if (s_l<=0) return 0.0; else return s_l*s_l;
  }

  //! gas phase relative permeability
  typename Traits::RangeFieldType 
  kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType s_g) const
  {
    if (s_g<=0) return 0.0; else return s_g*s_g;
  }

  //! liquid phase viscosity
  typename Traits::RangeFieldType 
  mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_l) const
  {
    return 1E-3;
  }

  //! gas phase viscosity
  typename Traits::RangeFieldType 
  mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_g) const
  {
    return 1E-5;
  }
	  
  //! absolute permeability (scalar!)
  typename Traits::RangeFieldType 
  k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1E-10;
  }

  //! gravity vector
  const typename Traits::RangeType& gravity () const
  {
    return gvector;
  }

  //! liquid phase molar density
  typename Traits::RangeFieldType 
  nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_l) const
  {
    return 1000.0;
  }

  //! gas phase molar density
  typename Traits::RangeFieldType 
  nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
        typename Traits::RangeFieldType p_g) const
  {
    return 1.0;
  }

  //! liquid phase mass density
  typename Traits::RangeFieldType 
  rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
         typename Traits::RangeFieldType p_l) const
  {
    return 1000.0;
  }

  //! gas phase mass density
  typename Traits::RangeFieldType 
  rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
         typename Traits::RangeFieldType p_g) const
  {
    return 1.0;
  }
	  
  //! liquid phase boundary condition type
  int
  bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 1;
  }

  //! gas phase boundary condition type
  int
  bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 1;
  }

  //! liquid phase Dirichlet boundary condition
  typename Traits::RangeFieldType 
  g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 1.0;
  }

  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType 
  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 1.0;
  }

  //! liquid phase Neumann boundary condition
  typename Traits::RangeFieldType 
  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  //! gas phase Neumann boundary condition
  typename Traits::RangeFieldType 
  j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  //! liquid phase source term
  typename Traits::RangeFieldType 
  q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }
  
  //! gas phase source term
  typename Traits::RangeFieldType 
  q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

private:
  typename Traits::RangeType gvector;
};

//==============================================================================
// saturation output
//==============================================================================

template<typename  T, typename PL, typename PG>
class S_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_l<T,PL,PG> >
{
  const T& t;
  const PL pl;
  const PG pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
    typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
    typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_l<T,PL,PG> > BaseT;

  S_l (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = t.s_l(e,x,pg_value-pl_value);
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }
};

template<typename  T, typename PL, typename PG>
class S_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_g<T,PL,PG> >
{
  const T& t;
  const PL pl;
  const PG pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
    typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
    typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_g<T,PL,PG> > BaseT;

  S_g (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = 1.0-t.s_l(e,x,pg_value-pl_value);
  }
  
  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }
};


//==============================================================================
// driver
//==============================================================================

template<class GV> 
void test (const GV& gv)
{
  // some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType::cube);
  
  // make grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::NoConstraints,Dune::PDELab::ISTLVectorBackend<2>,
    Dune::PDELab::SimpleGridFunctionStaticSize> GFS;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,
    Dune::PDELab::GridFunctionSpaceBlockwiseMapper> TPGFS;
  watch.reset();
  GFS gfs(gv,fem);
  TPGFS tpgfs(gfs);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // initial value function
  typedef P_l<GV,RF> P_lType;
  P_lType p_l_initial(gv);
  typedef P_g<GV,RF> P_gType;
  P_gType p_g_initial(gv);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,P_gType> PType;
  PType p_initial(p_l_initial,p_g_initial);

  // make vector for old time step and initialize
  typedef typename TPGFS::template VectorContainer<RF>::Type V;
  V pold(tpgfs);
  Dune::PDELab::interpolate(p_initial,tpgfs,pold);

  // make vector for new time step and initialize
  V pnew(tpgfs);
  pnew = pold;

  // make local operator
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp;
  typedef Dune::PDELab::TwoPhaseTwoPointFluxOperator<TP,V> LOP;
  LOP lop(tp,pold);
  lop.set_time(0.0);
  lop.set_timestep(0.001);

  // make grid operator space
  typedef typename TPGFS::template ConstraintsContainer<RF>::Type ET;
  ET et;
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<2,2> MB;
  typedef Dune::PDELab::GridOperatorSpace<TPGFS,TPGFS,LOP,ET,ET,MB> TPGOS;
  TPGOS tpgos(tpgfs,et,tpgfs,et,lop);

  // represent operator as a matrix
  typedef typename TPGOS::template MatrixContainer<RF>::Type M;
  M m(tpgos);
  m = 0.0;
  watch.reset();
  tpgos.jacobian(pnew,m);
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // compute residual
  V r(tpgfs,0.0);
  tpgos.residual(pnew,r);

  // make subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,0> P_lSUB;
  P_lSUB p_lsub(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,1> P_gSUB;
  P_gSUB p_gsub(tpgfs);
 
  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<P_lSUB,V> P_lDGF;
  P_lDGF p_ldgf(p_lsub,pnew);
  typedef Dune::PDELab::DiscreteGridFunction<P_gSUB,V> P_gDGF;
  P_gDGF p_gdgf(p_gsub,pnew);
  typedef S_l<TP,P_lDGF,P_gDGF> S_lDGF; 
  S_lDGF s_ldgf(tp,p_ldgf,p_gdgf);
  typedef S_g<TP,P_lDGF,P_gDGF> S_gDGF; 
  S_gDGF s_gdgf(tp,p_ldgf,p_gdgf);

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_lDGF>(p_ldgf,"p_l"));
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<P_gDGF>(p_gdgf,"p_g"));
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_gDGF>(s_gdgf,"s_g"));
  vtkwriter.write("twophase",Dune::VTKOptions::ascii);
}

//==============================================================================
// grid setup
//==============================================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(3);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      //grid.globalRefine(1);
      
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
