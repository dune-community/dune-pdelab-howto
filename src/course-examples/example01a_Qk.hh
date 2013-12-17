template<int k, class GV>
void example01a_Qk (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,k> FEM;
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  // <<<3>>> Make grid operator
  typedef Example01aLocalOperator LOP;
  LOP lop;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  // Structured 2D grid, Q1 finite elements => 9-point stencil / Q2 => 25
  MBE mbe(k == 1 ? 9 : 25);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,gfs,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  typename GO::Traits::Jacobian jac(go);
  std::cout << jac.patternStatistics() << std::endl;

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
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
  vtkwriter.write("example01_Q1",Dune::VTK::appendedraw);
}
