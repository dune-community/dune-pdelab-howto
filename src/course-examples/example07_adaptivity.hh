template<class Grid, class GV>
void adaptivity (Grid& grid, const GV& gv, int startLevel, int maxLevel)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,Coord,Real,1> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  BCTypeParam bctype; // boundary condition type
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc );               // assemble constraints

  // <<<3>>> Make grid operator
  typedef Example02LocalOperator<BCTypeParam> LOP;       // operator including boundary
  LOP lop(bctype);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(7);
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  typename GO::Traits::Jacobian jac(go);
  std::cout << jac.patternStatistics() << std::endl;

  // <<<4>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef BCExtension<GV,Real> G;                        // boundary value + extension
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,u);                    // interpolate coefficient vector

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> Assemble linear problem.
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,ls,u,1e-10);


  // <<<7>>> Preparation: Define types for the computation of the error estimate eta.
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> P0FEM;
  P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,VBE> P0GFS;
  typedef Dune::PDELab::ExampleErrorEstimator ESTLOP;
  typedef Dune::PDELab::EmptyTransformation NoTrafo;
  typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,Real>::Type U0;

  for (int i = 0; i <= maxLevel - startLevel; i++)
  {
    std::stringstream s;
    s << i;
    std::string iter;
    s >> iter;
    std::cout << "Iteration: " << iter << "\thighest level in grid: " << grid.maxLevel() << std::endl;
    std::cout << "constrained dofs=" << cc.size()
            << " of " << gfs.globalSize() << std::endl;

    // <<<8>>> Solve linear problem.
    slp.apply();

    // <<<9>>> graphical output
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
    vtkwriter.write("adaptivity_"+iter,Dune::VTK::appendedraw);

    // <<<10>>> compute estimated error eta
    P0GFS p0gfs(gv,p0fem);
    ESTLOP estlop;
    typedef Dune::PDELab::GridOperator<GFS,P0GFS,ESTLOP,MBE,Real,Real,Real,NoTrafo,NoTrafo> ESTGO;
    ESTGO estgo(gfs,p0gfs,estlop,mbe);
    U0 eta(p0gfs,0.0);
    estgo.residual(u,eta);

    for (unsigned int i=0; i<eta.flatsize(); i++)
      eta.base()[i] = sqrt(eta.base()[i]);

    // Use eta to refine the grid following two different strategies based
    // (1) element fraction
    // (2) error fraction

    double alpha(0.4);       // refinement fraction
    double eta_alpha(0);     // refinement threshold
    double beta(0.0);        // coarsening fraction
    double eta_beta(0);      // coarsening threshold
    int verbose = 2;

    // <<<10>>> Adapt the grid locally...
    // with strategy 1:
    element_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );
    // or, alternatively, with strategy 2:
    //error_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );

    mark_grid( grid, eta, eta_alpha, 0.0 );
    adapt_grid( grid, gfs, u, 2 );

    Dune::PDELab::constraints(bctype,gfs,cc);
    Dune::PDELab::interpolate(g,gfs,u);
  }
}
