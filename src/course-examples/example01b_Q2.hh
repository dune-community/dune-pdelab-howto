template<class GV>
void example01b_Q2 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,2> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  // <<<3>>> Make DOF vector
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0); // initial value

  // <<<4>>> Make grid operator
  typedef Example01bLocalOperator LOP;                                     // <= NEW
  LOP lop(4);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(25);

  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,gfs,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  typename GO::Traits::Jacobian jac(go);
  std::cout << jac.patternStatistics() << std::endl;

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,U> newton(go,u,ls);                         // <= NEW
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(10);
  newton.apply();

  // <<<7>>> graphical output
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
  vtkwriter.write("example01b_Q2",Dune::VTK::appendedraw);       // <= NEW
}
