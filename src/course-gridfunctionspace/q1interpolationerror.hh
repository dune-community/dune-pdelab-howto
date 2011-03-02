#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include"l2interpolationerror.hh"

template<typename GV>
void q1interpolationerror (const GV& gv)
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  Q1LocalFiniteElementMap<D,R> fem;   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<GV,
	Q1LocalFiniteElementMap<D,R> > GFS;    
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename Dune::PDELab::BackendVectorSelector<GFS,R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // make analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // make x interpolate u

  std::cout.precision(8);
  std::cout << "interpolation error: " 
			<< std::setw(8) << gv.size(0) << " elements " 
			<< std::scientific << l2interpolationerror(u,gfs,x,4) << std::endl;
}
