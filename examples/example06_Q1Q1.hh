template<class GV>
void example06_Q1Q1 (const GV& gv, double dtstart, double dtmax, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> initialize time variable
  Real time = 0.0;

  // <<<2a>>> Make grid function space for the system
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM0;
  FEM0 fem0;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
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

  // <<<2b>>> Make subspaces to extract components for visualization
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,0> U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace<GFS,1> U1SUB;
  U1SUB u1sub(gfs);

  // <<<3>>> Constraints on function space
  typedef BCType<GV> U0BC;
  U0BC u0bc(gv);
  typedef Dune::PDELab::PowerGridFunction<U0BC,2> UBC;
  UBC ubc(u0bc);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(ubc,gfs,cc);
  if (gv.comm().rank()<0) std::cout << "/" << gv.comm().rank() << "/ " << "constrained dofs=" << cc.size() 
				    << " of " << gfs.globalSize() << std::endl;

  // <<<4>>> Make FE function with initial value
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V xold(gfs,0.0);
  typedef U0Initial<GV,Real> U0InitialType;
  U0InitialType u0initial(gv);
  typedef U1Initial<GV,Real> U1InitialType;
  U1InitialType u1initial(gv);
  typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType;
  UInitialType uinitial(u0initial,u1initial);
  Dune::PDELab::interpolate(uinitial,gfs,xold);

  // <<<4>>> Make instationary grid operator space
  Real d_0 = 0.00028;
  Real d_1 = 0.005;
  Real lambda = 1.0;
  Real sigma = 1.0;
  Real kappa = -0.05;
  Real tau = 0.1;
  typedef Example05LocalOperator LOP; 
  LOP lop(d_0,d_1,lambda,sigma,kappa,4);
  typedef Example05TimeLocalOperator TLOP; 
  TLOP tlop(tau,4);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<2,2> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,TLOP,CC,CC,MBE> IGOS;
  IGOS igos(gfs,cc,gfs,cc,lop,tlop);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(gfs,cc,5000,5,1);

  // <<<6>>> Solver for non-linear problem per stage
  typedef Dune::PDELab::Newton<IGOS,LS,V> PDESOLVER;
  PDESOLVER pdesolver(igos,ls);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setReduction(1e-10);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);

  // <<<7>>> time-stepper
  Dune::PDELab::FractionalStepParameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,V,V> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example06_Q1Q1");
  {
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
    U0DGF u0dgf(u0sub,xold);
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
    U1DGF u1dgf(u1sub,xold);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<9>>> time loop
  V xnew(gfs,0.0);
  xnew = xold;
  double dt = dtstart;
  while (time<tend-1e-8)
    {
      // do time step
      osm.apply(time,dt,xold,xnew);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
      U0DGF u0dgf(u0sub,xnew);
      typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
      U1DGF u1dgf(u1sub,xnew);
      Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
      fn.increment();

      xold = xnew;
      time += dt;
      if (dt<dtmax-1e-8)
        dt = std::min(dt*1.2,dtmax);
    }
}
