#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>

template<typename GV> void q1GridFunctionSpace (const GV& gv) {
  typedef typename GV::Grid::ctype D;// domain type
  typedef double R;                  // range type

  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,D,R,1> FEM;
  FEM fem(gv);                           // maps entity to finite element

  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; /*@\label{q1gfs:GFS}@*/
  GFS gfs(gv,fem);                   // make grid function space

  typedef typename Dune::PDELab::BackendVectorSelector<GFS,R>::Type X;
  X x(gfs,0.0);                      // make coefficient vector

  typedef Dune::PDELab::DiscreteGridFunction<GFS,X> DGF;
  DGF dgf(gfs,x);                    // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);  // plot result
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"q1"));
  vtkwriter.write("q1gridfunctionspace",Dune::VTK::ascii);
}
