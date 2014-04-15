// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Simple Maxwell problem in 3D
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

#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
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
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/localoperator/maxwelldg.hh>

#include"../utility/gridexamples.hh"

//==============================================================================
// Parameter class for Maxwell Problem
//==============================================================================

template<typename GV, typename RF>
class MaxwellModelProblem
{
public:
  typedef Dune::PDELab::MaxwellParameterTraits<GV,RF> Traits;

  MaxwellModelProblem ()
    : pi(3.141592653589793238462643), time(0.0)
  {
  }

  //! permittivity
  typename Traits::RangeFieldType
  eps (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  //! permeability
  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  //! permeability
  typename Traits::RangeFieldType
  sigma (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition value
  typename Traits::StateType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::StateType u(0.0);
    // u[0] = -s[0];
    // u[1] = -s[1];
    // u[2] = -s[2];
    // u[3] = -s[3];
    // u[4] = -s[4];
    // u[5] = -s[5];
    return u;
  }

  //! right hand side
  typename Traits::StateType
  j (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& xin) const
  {
    typedef typename Traits::DomainFieldType DF;
    typename Traits::DomainType xglobal = e.geometry().global(xin);
    DF x=xglobal[0], y=xglobal[1], z=xglobal[2];
    typename Traits::StateType u(0.0);
    if (x>0.4 && x<0.6 && y>0.4 && y<0.6 && z>0.4 && z<0.6)
      {
        RF p1 = sin(pi*(x-0.4)/0.2)*sin(pi*(x-0.4)/0.2);
        RF p1prime = 2.0*sin(pi*(x-0.4)/0.2)*cos(pi*(x-0.4)/0.2)*pi/0.2;
        RF q2 = sin(pi*(y-0.4)/0.2)*sin(pi*(y-0.4)/0.2);
        RF q2prime = 2.0*sin(pi*(y-0.4)/0.2)*cos(pi*(y-0.4)/0.2)*pi/0.2;
        RF r3 = sin(pi*(z-0.4)/0.2)*sin(pi*(z-0.4)/0.2);
        RF r3prime = 2.0*sin(pi*(z-0.4)/0.2)*cos(pi*(z-0.4)/0.2)*pi/0.2;
        u[0] = -p1*q2prime*r3prime;
        u[1] = -0.5*p1prime*q2*r3prime;
        u[2] = -0.5*p1prime*q2prime*r3;
      }
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

//===============================================================
// driver
//===============================================================

// example using explicit time-stepping
template<class GV, class FEMDG, int degree>
void explicit_scheme (const GV& gv, const FEMDG& femdg, double Tend, double timestep, std::string name, int modulo)
{
  std::cout << "using degree " << degree << std::endl;
  // <<<1>>> Choose domain and range field type
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend
    <Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE> GFSDG;
  GFSDG gfsdg(gv,femdg);
  typedef Dune::PDELab::PowerGridFunctionSpace
    <GFSDG,dim*2,Dune::PDELab::ISTLVectorBackend<> > GFS;
  GFS gfs(gfsdg);
  gfs.update(); // required here for initialization of the gfs
  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  std::cout << "degrees of freedom: " << gfs.globalSize() << std::endl;

  // <<<2b>>> define problem parameters
  typedef MaxwellModelProblem<GV,Real> Param;
  Param param;

  // <<<3>>> Make grid operator space
  typedef Dune::PDELab::DGMaxwellSpatialOperator<Param,FEMDG> LOP;
  LOP lop(param);
  typedef Dune::PDELab::DGMaxwellTemporalOperator<Param,FEMDG> TLOP;
  TLOP tlop(param);
  Dune::PDELab::ExplicitEulerParameter<Real> method1;
  Dune::PDELab::HeunParameter<Real> method2;
  Dune::PDELab::Shu3Parameter<Real> method3;
  Dune::PDELab::RK4Parameter<Real> method4;
  Dune::PDELab::TimeSteppingParameterInterface<Real> *method;
  if (degree==0) {method=&method1; std::cout << "setting explicit Euler" << std::endl;}
  if (degree==1) {method=&method2; std::cout << "setting Heun" << std::endl;}
  if (degree==2) {method=&method3; std::cout << "setting Shu 3" << std::endl;}
  if (degree==3) {method=&method4; std::cout << "setting RK4" << std::endl;}

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  typedef Dune::PDELab::GridOperator
    <GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO0;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);

  typedef Dune::PDELab::GridOperator
    <GFS,GFS,TLOP,MBE,Real,Real,Real,C,C> GO1;
  GO1 go1(gfs,cg,gfs,cg,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
  IGO igo(go0,go1);
  igo.setMethod(*method);

  // <<<4>>> set initial values
  typedef typename IGO::Traits::Domain V;
  V xold(gfs,0.0);
  Dune::PDELab::MaxwellInitialValueAdapter<Param> u0(gv,param);
  Dune::PDELab::interpolate(u0,gfs,xold);

  // <<<5>>> Make a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  //LS ls(10000,1);
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<6>>> time-stepper
  typedef Dune::PDELab::SimpleTimeController<Real> TC;
  TC tc;
  Dune::PDELab::ExplicitOneStepMethod<Real,IGO,LS,V,V,TC> osm(*method,igo,ls,tc);
  osm.setVerbosityLevel(2);

  // <<<7>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn(name);
  using namespace Dune::TypeTree;
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<0> > U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<1> > U1SUB;
  U1SUB u1sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<2> > U2SUB;
  U2SUB u2sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<3> > U3SUB;
  U3SUB u3sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<4> > U4SUB;
  U4SUB u4sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,TreePath<5> > U5SUB;
  U5SUB u5sub(gfs);
  {
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
    U0DGF u0dgf(u0sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
    U1DGF u1dgf(u1sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U2SUB,V> U2DGF;
    U2DGF u2dgf(u2sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U3SUB,V> U3DGF;
    U3DGF u3dgf(u3sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U4SUB,V> U4DGF;
    U4DGF u4dgf(u4sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U5SUB,V> U5DGF;
    U5DGF u5dgf(u5sub,xold);
    int refinement = std::max(degree-1,0);
    if (degree>=2) refinement+=1;
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"D_x"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"D_y"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U2DGF>(u2dgf,"D_z"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U3DGF>(u3dgf,"B_x"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U4DGF>(u4dgf,"B_y"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U5DGF>(u5dgf,"B_z"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTK::appendedraw);
    fn.increment();
  }

  // <<<8>>> time loop
  int counter=0;
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
          typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
          U0DGF u0dgf(u0sub,x);
          typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
          U1DGF u1dgf(u1sub,x);
          typedef Dune::PDELab::DiscreteGridFunction<U2SUB,V> U2DGF;
          U2DGF u2dgf(u2sub,x);
          typedef Dune::PDELab::DiscreteGridFunction<U3SUB,V> U3DGF;
          U3DGF u3dgf(u3sub,x);
          typedef Dune::PDELab::DiscreteGridFunction<U4SUB,V> U4DGF;
          U4DGF u4dgf(u4sub,x);
          typedef Dune::PDELab::DiscreteGridFunction<U5SUB,V> U5DGF;
          U5DGF u5dgf(u5sub,x);
          int refinement = std::max(degree-1,0);
          if (degree>=2) refinement+=1;
          Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"D_x"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"D_y"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U2DGF>(u2dgf,"D_z"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U3DGF>(u3dgf,"B_x"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U4DGF>(u4dgf,"B_y"));
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U5DGF>(u5dgf,"B_z"));
          vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTK::appendedraw);
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

    if (argc!=7)
      {
        if(helper.rank()==0)
          {
            std::cout << "usage: " << argv[0] << " <end time> <time step> <grid file> <refinement> <degree> <modulo>" << std::endl;
            std::cout << "         <grid file> = 'yaspgrid' || <a gmsh file>"  << std::endl;
            std::cout << "         <refinement> = #cell per dir in yaspgrid, #refinements in UG" << std::endl;
            std::cout << "         <modulo> = write vtk file every modulo'th time step" << std::endl;
          }
        return 1;
      }

    double Tend;
    sscanf(argv[1],"%lg",&Tend);
    double timestep;
    sscanf(argv[2],"%lg",&timestep);
    std::string grid_file(argv[3]);
    int max_level; sscanf(argv[4],"%d",&max_level);
    int p; sscanf(argv[5],"%d",&p);
    int modulo; sscanf(argv[6],"%d",&modulo);

    // parallel overlapping yaspgrid version
    if (grid_file=="yaspgrid")
      {
        const int dim = 3;
        Dune::FieldVector<double,dim> L(1.0);
        Dune::array<int,dim> N(Dune::fill_array<int,dim>(max_level));
        std::bitset<dim> periodic(false);
        int overlap=1;
        Dune::YaspGrid<dim> grid(helper.getCommunicator(),L,N,periodic,overlap);
        typedef Dune::YaspGrid<dim>::LeafGridView GV;
        const GV& gv=grid.leafGridView();
        if (p==0)
          {
            const int degree=0;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_n" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==1)
          {
            const int degree=1;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_n" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==2)
          {
            const int degree=2;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_n" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        // if (p==3)
        //   {
        //     const int degree=3;
        //     typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
        //     FEM fem;
        //     std::stringstream fullname;
        //     fullname << grid_file << "_l" << max_level << "_k" << p;
        //     explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
        //   }
        return 0;
      }

#if HAVE_UG
    if (true)
      {
        // make uggrid
        const int dim=3;
        typedef Dune::UGGrid<dim> GridType;
        typedef std::vector<int> GmshIndexMap;
        GmshIndexMap boundary_index_map;
        GmshIndexMap element_index_map;
        Dune::GmshReader<GridType> gmsh_reader;
        Dune::shared_ptr<GridType>
          gridp(gmsh_reader.read(grid_file,boundary_index_map,
                                 element_index_map,true,false));
        for (int i=0; i<max_level; i++) gridp->globalRefine(1);
        typedef GridType::LeafGridView GV;
        const GV& gv=gridp->leafGridView();
        if (p==0)
          {
            const int degree=0;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==1)
          {
            const int degree=1;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==2)
          {
            const int degree=2;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        // if (p==3)
        //   {
        //     const int degree=3;
        //     typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
        //     FEM fem;
        //     std::stringstream fullname;
        //     fullname << grid_file << "_l" << max_level << "_k" << p;
        //     explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
        //   }
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
