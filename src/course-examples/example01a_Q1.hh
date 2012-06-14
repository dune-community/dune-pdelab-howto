template<class GV>
void example01a_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  // <<<3>>> Make grid operator
  typedef Example01aLocalOperator LOP; 
  LOP lop;
  typedef VBE::MatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,gfs,lop);

  // <<<4>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<5>>> solve linear problem
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0); // initial value
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,u,ls,1e-10);
  slp.apply();

  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("example01_Q1",Dune::VTK::appendedraw);
}
