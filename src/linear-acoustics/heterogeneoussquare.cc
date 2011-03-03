// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Acoustic wave propagation in a simple 2D cavity

    Example 4 from "L. Krivodonova, J. Xin, J.-F. Remacle, N. Chevaugeon, J.E. Flaherty:
    Shock detection and limiting with discontinuous Galerkin methods for hyperbolic
    conservation laws. Applied Numerical Mathematics, 48, 323-338, 2004.
*/
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

#include<dune/pdelab/finiteelementmap/opbfem.hh>
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
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/localoperator/linearacousticsdg.hh>

#include"../utility/gridexamples.hh"

//==============================================================================
// Parameter class for the linear acoustics problem
//==============================================================================

template<typename GV, typename RF>
class Krivodonova4Problem
{
public:
  typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,RF> Traits;

  Krivodonova4Problem ()
    : pi(3.141592653589793238462643), time(0.0)
  {
  }

  //! speed of sound
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 340.0;
  }

  //! Dirichlet boundary condition value
  typename Traits::StateType 
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1e-6) 
      {
        typename Traits::StateType u(0.0);
        u[1] = 1.224*(1+0.5*sin(2*pi*1500.0*time));
        return u;
      }
    if (xglobal[0]>1.0-1e-6) 
      {
        typename Traits::StateType u(0.0);
        return u;
      }
    typename Traits::StateType u(0.0);
    u[2] = 0.0;
    return u;
  }

  //! right hand side
  typename Traits::StateType 
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType 
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType u(0.0);
    return u;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  double pi;
  RF time;
};


template<typename GV, typename RF>
class RiemannProblem
{
public:
  typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,RF> Traits;

  RiemannProblem ()
    : time(0.0)
  {
  }

  //! speed of sound
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if (xglobal[1]<0.5) return 1.0;
    if (xglobal[0]>0.5) return 2.0;
    return 0.5;
  }

  //! Dirichlet boundary condition value
  typename Traits::StateType 
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::StateType u(0.0);
    u[0] = s[0];
    u[1] = -s[1];
    u[2] = -s[2];
    return u;
  }

  //! right hand side
  typename Traits::StateType 
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType 
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::StateType u(0.0);
    if (xglobal[0]>0.45 && xglobal[0]<0.55 && xglobal[1]>0.3 && xglobal[1]<0.4)
      {
        u[0] = 1;
        u[1] = 0;
        u[2] = 0;
      }
    return u;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  RF time;
};


//===============================================================
// Some variants to solve the nonlinear diffusion problem
//===============================================================


// example using explicit time-stepping
template<class GV>
void explicit_scheme (const GV& gv, double Tend, double timestep)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  const int degree=2;
  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
  typedef Dune::PDELab::OPBLocalFiniteElementMap<typename GV::Grid::ctype,Real,degree,dim,Dune::GeometryType::cube> FEMDG;
  FEMDG femdg;
  typedef Dune::PDELab::NoConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE,Dune::PDELab::SimpleGridFunctionStaticSize> GFSDG;
  GFSDG gfsdg(gv,femdg);
  typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFMapper;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFSDG,dim+1,GFMapper> GFS;
  GFS gfs(gfsdg);
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;

  // <<<2b>>> define problem parameters
  typedef RiemannProblem<GV,Real> Param;
  Param param;

  // <<<4>>> set initial values
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V xold(gfs,0.0);
  Dune::PDELab::LinearAcousticsInitialValueAdapter<Param> u0(gv,param);
  Dune::PDELab::interpolate(u0,gfs,xold);

  // <<<5>>> Make grid operator space
  typedef Dune::PDELab::DGLinearAcousticsSpatialOperator<Param,FEMDG> LOP; 
  LOP lop(param);
  typedef Dune::PDELab::DGLinearAcousticsTemporalOperator<Param,FEMDG> TLOP; 
  TLOP tlop(param);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<blocksize,blocksize> MBE;
  //Dune::PDELab::ExplicitEulerParameter<Real> method;
  Dune::PDELab::HeunParameter<Real> method;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,TLOP,C,C,MBE> IGOS;
  IGOS igos(method,gfs,cg,gfs,cg,lop,tlop);

  // <<<6>>> Make a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  //LS ls(10000,1);
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<8>>> time-stepper
  typedef Dune::PDELab::CFLTimeController<Real,IGOS> TC;
  TC tc(0.999,igos);
  Dune::PDELab::ExplicitOneStepMethod<Real,IGOS,LS,V,V,TC> osm(method,igos,ls,tc);
  osm.setVerbosityLevel(2);

  // <<<10>>> graphics for initial guess
  std::stringstream fullname;
  fullname << "riemann" << "_dim" << dim << "_k" << degree;
  Dune::PDELab::FilenameHelper fn(fullname.str());
  int counter=0, modulo=20;
  {
    typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    int refinement = std::max(degree-1,0);
    if (degree>=2) refinement+=2;
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"u"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  V x(gfs,0.0);
  while (time < Tend)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // graphics
      counter++;
      if (counter%modulo==0)
        {
          typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
          DGF xdgf(gfs,x);
          int refinement = std::max(degree-1,0);
          if (degree>=2) refinement+=2;
          Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"u"));
          vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTKOptions::binaryappended);
          fn.increment();
        }

      xold = x;
      time += dt;
    }
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

	if (argc!=5)
	  {
		if(helper.rank()==0)
		  std::cout << "usage: ./transporttest <end time> <time step> <elements in y> <overlap>" << std::endl;
		return 1;
	  }

	double Tend;
	sscanf(argv[1],"%lg",&Tend);

	double timestep;
	sscanf(argv[2],"%lg",&timestep);

	int n;
	sscanf(argv[3],"%d",&n);

	int o;
	sscanf(argv[4],"%d",&o);

    // parallel overlapping version
    if (true)
      {
        Dune::FieldVector<double,2> L; L[0]=1.0; L[1]=1.0;
        Dune::FieldVector<int,2> N; N[0]=n; N[1]=n;
        Dune::FieldVector<bool,2> periodic(false);
        int overlap=o;
        Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafView();
        explicit_scheme(gv,Tend,timestep);
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
