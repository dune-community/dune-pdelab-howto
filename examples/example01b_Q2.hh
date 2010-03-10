template<class GV>
void example01b_Q2 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2a>>> Make grid function space
  typedef Dune::PDELab::Q22DLocalFiniteElementMap<Coord,Real> FEM;
  FEM fem;
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  // <<<2b>>> Constraints (are empty here)
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;

  // <<<3>>> Make grid operator space
  typedef Example01bLocalOperator LOP; 
  LOP lop(4);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,gfs,lop);

  // <<<4>>> Select a linear solver backend
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  LS ls(true);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);
#endif

  // <<<5>>> solve nonlinear problem
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U x(gfs,2.0); // initial value
  Dune::PDELab::Newton<GOS,LS,U> newton(gos,x,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(10);
  newton.apply();

  // <<<6>>> graphical output
  {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF xdgf(gfs,x);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"solution"));
    vtkwriter.write("example01b_Q2",Dune::VTKOptions::binaryappended);
  }
}
