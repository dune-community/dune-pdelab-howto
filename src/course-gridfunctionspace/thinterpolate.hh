#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>

template<class GV> 
void thinterpolate (const GV& gv)
{
  // types
  typedef typename GV::Grid::ctype D;
  typedef double R;
  const int dim = GV::dimension;

  typedef Dune::PDELab::ISTLVectorBackend<> VBE;

  //  make Q_1 grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,D,R,1> Q1FEM;
  Q1FEM q1fem(gv);                    // Q1 finite elements
  typedef Dune::PDELab::GridFunctionSpace<GV,Q1FEM,VBE> Q1GFS;
  Q1GFS q1gfs(gv,q1fem);              // Q1 space
  q1gfs.name("pressure");
  
  // make Q_2 finite element map
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,D,R,2> Q2FEM;
  Q2FEM q2fem(gv);                    // Q2 finite elements

  // make velocity grid function space
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,Q2FEM,dim,VBE,VBE> VGFS;
  VGFS vgfs(gv,q2fem);                // velocity space
  vgfs.name("velocity");

  // make Taylor-Hood grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<VBE,
    Dune::PDELab::LexicographicOrderingTag,
    VGFS,Q1GFS> THGFS;              
  THGFS thgfs(vgfs,q1gfs);            // Taylor-Hood space

  // make coefficent vector
  typedef typename Dune::PDELab::BackendVectorSelector<THGFS,R>::Type X;
  X x(thgfs,0.0);                     // one x for all dofs !

  // interpolate from analytic function
  typedef U<GV,R> Pressure;
  Pressure p(gv);                     // pressure component
  typedef V<GV,R> Velocity;
  Velocity v(gv);                     // velocity component
  typedef Dune::PDELab::CompositeGridFunction<Velocity,Pressure> THF;
  THF thf(v,p);
  Dune::PDELab::interpolate(thf,thgfs,x);

  // output grid functions with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,thgfs,x);
  vtkwriter.write("thinterpolate",Dune::VTK::ascii);
}
