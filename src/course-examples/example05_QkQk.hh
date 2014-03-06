template<int k, class GV>
void example05_QkQk (const GV& gv, double dtstart, double dtmax, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  Real time = 0.0;

  // <<<2>>> Make grid function space for the system
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,k> FEM0;
  FEM0 fem0(gv);
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE0;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM0,CON,VBE0> GFS0;
  GFS0 gfs0(gv,fem0);

  typedef Dune::PDELab::ISTLVectorBackend
    <Dune::PDELab::ISTLParameters::static_blocking,2> VBE;               // block size 2
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS0,2,VBE,
    Dune::PDELab::EntityBlockedOrderingTag> GFS;
  GFS gfs(gfs0);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<0> > U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<1> > U1SUB;
  U1SUB u1sub(gfs);

  // <<<3>>> Make instationary grid operator
  Real d_0 = 0.00028;
  Real d_1 = 0.005;
  Real lambda = 1.0;
  Real sigma = 1.0;
  Real kappa = -0.05;
  Real tau = 0.1;
  typedef Example05LocalOperator LOP;
  LOP lop(d_0,d_1,lambda,sigma,kappa,2*k);
  typedef Example05TimeLocalOperator TLOP;
  TLOP tlop(tau,2*k);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(k == 1 ? 9 : 25);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO0;
  GO0 go0(gfs,gfs,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,Real,Real,Real,CC,CC> GO1;
  GO1 go1(gfs,gfs,tlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics (do not call this for IGO before osm.apply() was called!)
  typename GO0::Traits::Jacobian jac(go0);
  std::cout << jac.patternStatistics() << std::endl;

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
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);

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
  Dune::PDELab::Alexander2Parameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,U,U> osm(method,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
  std::stringstream basename;
  basename << "example05_Q" << k << "Q" << k;
  Dune::PDELab::FilenameHelper fn(basename.str());
  {
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,U> U0DGF;
    U0DGF u0dgf(u0sub,uold);
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,U> U1DGF;
    U1DGF u1dgf(u1sub,uold);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3*(k-1));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
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
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3*(k-1));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U0DGF>(u0dgf,"u0"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<U1DGF>(u1dgf,"u1"));
      vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
      fn.increment();

      uold = unew;
      time += dt;
      if (dt<dtmax-1e-8)
        dt = std::min(dt*1.1,dtmax);
    }
}
