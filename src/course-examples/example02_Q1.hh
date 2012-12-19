template<class GV>
void example02_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON; // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  BCTypeParam bctype; // boundary condition type
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc ); // assemble constraints
  std::cout << "constrained dofs=" << cc.size()
            << " of " << gfs.globalSize() << std::endl;

  // <<<3>>> Make grid operator
  typedef Example02LocalOperator<BCTypeParam> LOP;             // operator including boundary
  LOP lop( bctype );
  typedef Dune::PDELab::ISTLMatrixBackend MBE;

  typedef Dune::PDELab::GridOperator<
    GFS,GFS,        /* ansatz and test space */
    LOP,            /* local operator */
    MBE,            /* matrix backend */
    Real,Real,Real, /* field types for domain, range and jacobian */
    CC,CC           /* constraints transformation  for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop);

  // <<<4>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  typedef BCExtension<GV,Real> G;                               // boundary value + extension
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,u);                           // interpolate coefficient vector

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> assemble and solve linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,u,ls,1e-10);
  slp.apply();

  // <<<7>>> graphical output
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
  vtkwriter.write("example02_Q1",Dune::VTK::appendedraw);
}
