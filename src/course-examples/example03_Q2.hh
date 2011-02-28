template<class GV>
void example03_Q2 (const GV& gv, double dt, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  Real time = 0.0;                                              // make a time variable

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<Coord,Real> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  BCTypeParam bctype;
  bctype.setTime(time);                                              // b.c. depends on time now
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc );

  // <<<3>>> Make FE function with initial value / Dirichlet b.c.
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U uold(gfs,0.0);                                              // solution at t^n
  typedef BCExtension<GV,Real> G;                               // defines boundary condition,
  G g(gv);                                                      // extension and initial cond.
  g.setTime(time);                                              // b.c. depends on time now
  Dune::PDELab::interpolate(g,gfs,uold);

  // <<<4>>> Make instationary grid operator space
  typedef Example03LocalOperator<BCTypeParam> LOP; 
  LOP lop(bctype,4);                                           // local operator r
  typedef Example03TimeLocalOperator TLOP; 
  TLOP tlop(4);                                                 // local operator m
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,U,GFS,GFS,LOP,TLOP,CC,CC,MBE> IGOS;
  IGOS igos(gfs,cc,gfs,cc,lop,tlop);                            // new grid operator space

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  // <<<6>>> Solver for linear problem per stage
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGOS,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igos,ls,1e-10);

  // <<<7>>> time-stepper
  Dune::PDELab::Alexander2Parameter<Real> method;               // defines coefficients
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,U,U> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(2);                                     // time stepping scheme

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example03_Q2");              // append number to file name
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF udgf(gfs,uold);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
    fn.increment();                                             // increase file number
  }

  // <<<9>>> time loop
  U unew(gfs,0.0);                                              // solution to be computed
  while (time<tend-1e-8) {
      // do time step
      bctype.setTime(time+dt);                                       // compute constraints
      cc.clear();                                               // for this time step
      Dune::PDELab::constraints(bctype,gfs,cc);
      osm.apply(time,dt,uold,g,unew);                           // do one time step

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
      DGF udgf(gfs,unew);
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
      fn.increment();

      uold = unew;                                              // advance time step
      time += dt;
    }
}
