template<class GV>
void example05_Q1Q1 (const GV& gv, double dtstart, double dtmax, double tend) {
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  Real time = 0.0;

  // <<<2>>> Make grid function space for the system
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM0;
  FEM0 fem0;
  typedef Dune::PDELab::NoConstraints CON;                      // pure Neumann: no constraints
  typedef Dune::PDELab::ISTLVectorBackend<2> VBE;               // block size 2
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,CON,VBE> GFS0;
  GFS0 gfs0(gv,fem0);

  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM1;
  FEM1 fem1;                                                    // might use Q2 as well
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM1,CON,VBE> GFS1;
  GFS1 gfs1(gv,fem1);

  typedef Dune::PDELab::CompositeGridFunctionSpace<             // compose function space
  Dune::PDELab::GridFunctionSpaceBlockwiseMapper,GFS0,GFS1> GFS;// point block ordering
  GFS gfs(gfs0,gfs1);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  typedef Dune::PDELab::GridFunctionSubSpace<GFS,0> U0SUB;      // subspaces for later use
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,1> U1SUB;
  U1SUB u1sub(gfs);

  // <<<3>>> Make FE function with initial value
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U uold(gfs,0.0);
  typedef U0Initial<GV,Real> U0InitialType;
  U0InitialType u0initial(gv);
  typedef U1Initial<GV,Real> U1InitialType;
  U1InitialType u1initial(gv);
  typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType;
  UInitialType uinitial(u0initial,u1initial);
  Dune::PDELab::interpolate(uinitial,gfs,uold);

  // <<<4>>> Make instationary grid operator space
  Real d_0 = 0.00028, d_1 = 0.005, lambda = 1.0, sigma = 1.0, kappa = -0.05, tau = 0.1;
  typedef Example05LocalOperator LOP; 
  LOP lop(d_0,d_1,lambda,sigma,kappa,2);                        // spatial part
  typedef Example05TimeLocalOperator TLOP; 
  TLOP tlop(tau,2);                                             // temporal part
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<2,2> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,U,GFS,GFS,LOP,TLOP,CC,CC,MBE> IGOS;
  IGOS igos(gfs,gfs,lop,tlop);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

  // <<<6>>> Solver for non-linear problem per stage
  typedef Dune::PDELab::Newton<IGOS,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igos,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);

  // <<<7>>> time-stepper
  Dune::PDELab::Alexander2Parameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,U,U> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example05_Q1Q1");
  {
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
    U0DGF u0dgf(u0sub,uold);
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
    U1DGF u1dgf(u1sub,uold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<9>>> time loop
  U unew(gfs,0.0);
  unew = uold;
  double dt = dtstart;
  while (time<tend-1e-8)
    {
      // do time step
      osm.apply(time,dt,uold,unew);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
      U0DGF u0dgf(u0sub,unew);
      typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
      U1DGF u1dgf(u1sub,unew);
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
      fn.increment();

      uold = unew;
      time += dt;
      if (dt<dtmax-1e-8) dt = std::min(dt*1.1,dtmax);           // time step adaption
    }
}
