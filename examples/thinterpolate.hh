#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>

template<class GV> 
void thinterpolate (const GV& gv)
{
  // types
  typedef typename GV::Grid::ctype D;
  typedef double R;
  const int dim = GV::dimension;

  // make finite element maps
  typedef Dune::PDELab::Q12DLocalFiniteElementMap<D,R> Q12DFEM;
  Q12DFEM q12dfem;                    // Q1 finite elements
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<D,R> Q22DFEM;
  Q22DFEM q22dfem;                    // Q2 finite elements
  
  // make grid function spaces
  typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> Q1GFS;
  Q1GFS q1gfs(gv,q12dfem);            // Q1 space
  typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> Q2GFS;
  Q2GFS q2gfs(gv,q22dfem);            // Q2 space
  typedef Dune::PDELab::PowerGridFunctionSpace<Q2GFS,dim> VGFS;
  VGFS vgfs(q2gfs);                   // velocity space
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::GridFunctionSpaceLexicographicMapper,
	  VGFS,Q1GFS> THGFS;              
  THGFS thgfs(vgfs,q1gfs);            // Taylor-Hood space

  // make coefficent vector
  typedef typename THGFS::template VectorContainer<R>::Type X;
  X x(thgfs,0.0);                     // one x for all dofs !

  // make composite analytic function
  typedef U<GV,R> Pressure;
  Pressure p(gv);                     // pressure component
  typedef V<GV,R> Velocity;
  Velocity v(gv);                     // velocity component
  typedef Dune::PDELab::CompositeGridFunction<Velocity,Pressure> THF;
  THF thf(v,p);

  // do interpolation
  Dune::PDELab::interpolate(thf,thgfs,x);

  // select subspaces
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,0> VSUB;
  VSUB vsub(thgfs);                   // velocity subspace
  typedef Dune::PDELab::GridFunctionSubSpace<VSUB,0> V0SUB;
  V0SUB v0sub(vsub);                   // first velocity component
  typedef Dune::PDELab::GridFunctionSubSpace<THGFS,1> PSUB;
  PSUB psub(thgfs);                   // pressure subspace

  // make discrete function objects
  typedef Dune::PDELab::VectorDiscreteGridFunction<VSUB,X> VDGF;
  VDGF vdgf(vsub,x);
  typedef Dune::PDELab::DiscreteGridFunction<V0SUB,X> V0DGF;
  V0DGF v0dgf(v0sub,x);
  typedef Dune::PDELab::DiscreteGridFunction<PSUB,X> PDGF;
  PDGF pdgf(psub,x);

  // output grid functions with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(
	new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"velocity"));
  vtkwriter.addVertexData(
    new Dune::PDELab::VTKGridFunctionAdapter<V0DGF>(v0dgf,"velo 0"));
  vtkwriter.addVertexData(
    new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"pressure"));
  vtkwriter.write("thinterpolate",Dune::VTKOptions::ascii);
}
