template<class GV>
void example04 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem(Dune::GeometryType::cube);
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<3>>> Make grid operator space
  typedef Example04LocalOperator LOP; 
  LOP lop;
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::EmptyTransformation CC;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,gfs,lop);

  // <<<4>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  LS ls(5000,true);

  // <<<5>>> assemble and solve linear problem
  typedef typename GFS::template VectorContainer<Real>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,V> slp(gos,x,ls,1e-10);
  slp.apply();

  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF xdgf(gfs,x);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
  vtkwriter.write("example04",Dune::VTKOptions::binaryappended);
}
