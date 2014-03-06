template<class GV>
void example04 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));    // supply element type for P0
  typedef Dune::PDELab::NoConstraints CON;                      // we have no constraints!
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  // <<<3>>> Make grid operator
  typedef Example04LocalOperator LOP;                           // our new operator
  LOP lop;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5);
  typedef Dune::PDELab::EmptyTransformation CC;
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
  GO go(gfs,gfs,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  // typename GO::Traits::Jacobian jac(go);
  // std::cout << jac.patternStatistics() << std::endl;

  // <<<4>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  LS ls(100);

  // <<<5>>> assemble and solve linear problem
  typedef typename GO::Traits::Domain U;
  U u(gfs,0.0);
  Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> slp(go,ls,u,1e-10);
  slp.apply();

  // <<<6>>> graphical output
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,u);
  vtkwriter.write("example04",Dune::VTK::appendedraw);
}
