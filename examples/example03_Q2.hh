template<class GV>
void example03_Q2 (const GV& gv, double dt, double tend)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> initialize time variable
  Real time = 0.0;

  // <<<2a>>> Make grid function space
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<Coord,Real> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> Compute constraints on function space
  typedef BCType<GV> B;
  B b(gv);
  b.setTime(time); // depends on time now
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc);
  std::cout << "constrained dofs=" << cc.size() 
            << " of " << gfs.globalSize() << std::endl;

  // <<<3>>> Make FE function with initial value / Dirichlet b.c.
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V xold(gfs,0.0);
  typedef BCExtension<GV,Real> G;
  G g(gv);
  g.setTime(time);
  Dune::PDELab::interpolate(g,gfs,xold);

  // <<<4>>> Make instationary grid operator space
  typedef Example03LocalOperator<B> LOP; 
  LOP lop(b,4);
  typedef Example03TimeLocalOperator TLOP; 
  TLOP tlop(4);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::InstationaryGridOperatorSpace<Real,V,GFS,GFS,LOP,TLOP,CC,CC,MBE> IGOS;
  IGOS igos(gfs,cc,gfs,cc,lop,tlop);

  // <<<5>>> Select a linear solver backend
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls(false);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);
#endif

  // <<<6>>> Solver for linear problem per stage
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGOS,LS,V> PDESOLVER;
  PDESOLVER pdesolver(igos,ls,1e-10);

  // <<<7>>> time-stepper
  Dune::PDELab::Alexander2Parameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,IGOS,PDESOLVER,V,V> osm(method,igos,pdesolver);
  osm.setVerbosityLevel(2);

  // <<<8>>> graphics for initial guess
  Dune::PDELab::FilenameHelper fn("example03_Q2");
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
    fn.increment();
  }

  // <<<9>>> time loop
  V xnew(gfs,0.0);
  while (time<tend-1e-8)
    {
      // do time step
      b.setTime(time+dt);
      cc.clear();
      Dune::PDELab::constraints(b,gfs,cc);
      osm.apply(time,dt,xold,g,xnew);

      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
      DGF xdgf(gfs,xnew);
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
      vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
      fn.increment();

      xold = xnew;
      time += dt;
    }
}
