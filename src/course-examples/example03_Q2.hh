template<class GV>
void example03_Q2 (const GV& gv, double dt, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  Real time = 0.0;                                           // make a time variable

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,2> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  BCTypeParam bctype;
  bctype.setTime(time);                                      // b.c. depends on time now
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc );

  // <<<3>>> Make instationary grid operator
  typedef Example03LocalOperator<BCTypeParam> LOP; 
  LOP lop(bctype,4);                                         // local operator r
  typedef Example03TimeLocalOperator TLOP; 
  TLOP tlop(4);                                              // local operator m
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(25);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);                                          // new grid operator

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics (do not call this for IGO before osm.apply() was called!)
  typename GO0::Traits::Jacobian jac(go0);
  std::cout << jac.patternStatistics() << std::endl;

  // <<<3>>> Make FE function with initial value / Dirichlet b.c.
  typedef typename IGO::Traits::Domain U;
  U uold(gfs,0.0);                                           // solution at t^n
  typedef BCExtension<GV,Real> G;                            // defines boundary condition,
  G g(gv);                                                   // extension and initial cond.
  g.setTime(time);                                           // b.c. depends on time now
  Dune::PDELab::interpolate(g,gfs,uold);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  // <<<6>>> Solver for linear problem per stage
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igo,ls,1e-10);

  // <<<7>>> time-stepper
  Dune::PDELab::Alexander2Parameter<Real> method;          // defines coefficients
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,U,U> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(2);                                   // time stepping scheme

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example03_Q2");            // append number to file name
  {
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,uold);
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
    fn.increment();                                           // increase file number
  }

  // <<<9>>> time loop
  U unew(gfs,0.0);                                            // solution to be computed
  while (time<tend-1e-8) {
      // do time step
      bctype.setTime(time+dt);                                // compute constraints
      cc.clear();                                             // for this time step
      Dune::PDELab::constraints(bctype,gfs,cc);
      osm.apply(time,dt,uold,g,unew);                         // do one time step

      // graphics
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
      Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,unew);
      vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
      fn.increment();

      uold = unew;                                            // advance time step
      time += dt;
    }
}
