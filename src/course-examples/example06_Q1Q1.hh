template<class GV>
void example06_Q1Q1 (const GV& gv, double dtstart, double dtmax, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;
  Real time = 0.0;

  // <<<2>>> Make grid function space for the system
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM0;
  FEM0 fem0;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON; // new constraints class
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<2> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,CON,VBE> GFS0;
  GFS0 gfs0(gv,fem0,con);

  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM1;
  FEM1 fem1;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM1,CON,VBE> GFS1;
  GFS1 gfs1(gv,fem1,con);

  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::GridFunctionSpaceBlockwiseMapper,
	  GFS0,GFS1> GFS;              
  GFS gfs(gfs0,gfs1);

  typedef BCTypeParam<GV> U0_BCTypeParam;
  U0_BCTypeParam u0_bctype( gv );
  typedef Dune::PDELab::PowerConstraintsParameters<U0_BCTypeParam,2> U_BCTypeParam;
  U_BCTypeParam u_bctype( u0_bctype );
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;                                          // constraints needed due 
  Dune::PDELab::constraints( u_bctype, gfs, cc);  // to artificial boundaries

  typedef Dune::PDELab::GridFunctionSubSpace<GFS,0> U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,1> U1SUB;
  U1SUB u1sub(gfs);

  // <<<3>>> Make instationary grid operator space
  Real d_0 = 0.00028, d_1 = 0.005;
  Real lambda = 1.0, sigma = 1.0, kappa = -0.05, tau = 0.1;
  typedef Example05LocalOperator LOP; 
  LOP lop(d_0,d_1,lambda,sigma,kappa,2);
  typedef Example05TimeLocalOperator TLOP; 
  TLOP tlop(tau,2);
  typedef VBE::MatrixBackend MBE;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop);  
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);

  // <<<4>>> Make FE function with initial value
  typedef typename IGO::Traits::Domain U;
  U uold(gfs,0.0);
  typedef U0Initial<GV,Real> U0InitialType;
  U0InitialType u0initial(gv);
  typedef U1Initial<GV,Real> U1InitialType;
  U1InitialType u1initial(gv);
  typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType;
  UInitialType uinitial(u0initial,u1initial);
  Dune::PDELab::interpolate(uinitial,gfs,uold);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS; // select parallel backend !
  LS ls(gfs,cc,5000,5,1);

  // <<<6>>> Solver for non-linear problem per stage
  typedef Dune::PDELab::Newton<IGO,LS,U> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);

  // <<<7>>> time-stepper
  Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,U,U> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example06_Q1Q1");
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
      if (dt<dtmax-1e-8)
        dt = std::min(dt*1.1,dtmax);
    }
}
