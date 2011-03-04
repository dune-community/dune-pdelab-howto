template<class GV>
void example01a_Q2 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<Coord,Real> FEM;         // <= NEW
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  // <<<3>>> Make grid operator space
  typedef Example01aLocalOperator LOP; 
  LOP lop(4);                                                              // <= NEW
  typedef VBE::MatrixBackend MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,gfs,lop);

  // <<<4>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<5>>> solve linear problem
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0); // initial value
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  SLP slp(gos,u,ls,1e-10);
  slp.apply();

  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);                          // <= NEW
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("example01_Q2",Dune::VTKOptions::binaryappended);        // <= NEW
}
