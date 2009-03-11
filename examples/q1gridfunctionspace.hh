#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

template<typename GV>
void q1GridFunctionSpace (const GV& gv)
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  Q1LocalFiniteElementMap<D,R> fem;   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<GV,
	Q1LocalFiniteElementMap<D,R> > GFS;    
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector
  x[4] = 1.0;                         // set a component

  typedef Dune::PDELab::DiscreteGridFunction<GFS,X> DGF;
  DGF dgf(gfs,x);                     // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);  // plot result
  vtkwriter.addVertexData(new Dune::PDELab::
						  VTKGridFunctionAdapter<DGF>(dgf,"q1"));
  vtkwriter.write("q1gridfunctionspace",Dune::VTKOptions::ascii);
}
