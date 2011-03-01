template<class Grid, class GV>
void adaptivity (Grid& grid, const GV& gv, int startLevel, int maxLevel)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  BCTypeParam bctype; // boundary condition type
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints( bctype, gfs, cc );               // assemble constraints

  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U u(gfs,0.0);
  typedef BCExtension<GV,Real> G;                        // boundary value + extension
  G g(gv);
  Dune::PDELab::interpolate(g,gfs,u);                    // interpolate coefficient vector

  // <<<4>>> Make grid operator space
  typedef Example02LocalOperator<BCTypeParam> LOP;       // operator including boundary
  LOP lop(bctype);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
  
  // <<<6>>> assemble linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  SLP slp(gos,u,ls,1e-10);

  // <<<7>>> create GridAdaptor
  typedef Dune::PDELab::L2Projection<GFS,U> Proj;
  Proj proj;
  typedef Dune::PDELab::ResidualErrorEstimation<GFS,U> REE;
  REE ree(gfs);
  typedef Dune::PDELab::EstimationAdaptation<Grid,GFS,U,REE> EA;
  EA ea(grid,gfs,ree,0.2,0.4,1,maxLevel);
  typedef Dune::PDELab::GridAdaptor<Grid,GFS,U,EA,Proj> GRA;
  GRA gra(grid,gfs,ea,proj);

  for (int i = 0; i <= maxLevel - startLevel; i++)
  {
    std::stringstream s;
    s << i;
    std::string iter;
    s >> iter;
    std::cout << "Iteration: " << iter << "\thighest level in grid: " << grid.maxLevel() << std::endl;
  std::cout << "constrained dofs=" << cc.size() 
            << " of " << gfs.globalSize() << std::endl;
    
    // <<<8>>> solve linear problem
    slp.apply();

    // <<<9>>> graphical output
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF udgf(gfs,u);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write("adaptivity_"+iter,Dune::VTKOptions::binaryappended);
    
    // <<<10>>> adapt the grid
    gra.adapt(u);
    Dune::PDELab::constraints(bctype,gfs,cc);
    Dune::PDELab::interpolate(g,gfs,u);
  }
}
