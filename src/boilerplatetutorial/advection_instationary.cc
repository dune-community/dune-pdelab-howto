// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/l2.hh>

//***********************************************************************
//***********************************************************************
// diffusion problem with time dependent coefficients
//***********************************************************************
//***********************************************************************

template<typename GV, typename RF>
class GenericAdvectionProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericAdvectionProblem () : time(0.0)
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 0 : 0;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      v[i] = 0.0;
    v[0] = 1.0;
    v[1] = 1.0/3.0;
    b0 = 0.25;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    if (is.outerNormal(xlocal)*v<0)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    if (x[1] < v[1]*x[0]+b0)
      return 0.0;
    else
      return 1.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  typename Traits::PermTensorType I;
  typename Traits::RangeType v;
  typename Traits::RangeFieldType b0;
  RF time;
};


//***********************************************************************
//***********************************************************************
// a function that does the simulation on a given grid
//***********************************************************************
//***********************************************************************

template<typename GM, unsigned int degree, Dune::GeometryType::BasicType elemtype,
         Dune::PDELab::MeshType meshtype, Dune::SolverCategory::Category solvertype>
void do_simulation (double T, double dt, GM& grid, std::string basename)
{
  // define parameters
  typedef double NumberType;

  // make problem parameters
  typedef GenericAdvectionProblem<typename GM::LeafGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid.leafGridView(),problem);

  // make a finite element space
  typedef Dune::PDELab::DGPkSpace<GM,NumberType,degree,elemtype,solvertype> FS;
  FS fs(grid.leafGridView());

  // assemblers for finite element problem
  typedef Dune::PDELab::ConvectionDiffusionDG<Problem,typename FS::FEM> LOP;
  LOP lop(problem,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,2.0);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> SASS;
  SASS sass(fs,lop);
  typedef Dune::PDELab::L2 MLOP;
  MLOP mlop(2*degree);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,MLOP,solvertype> TASS;
  TASS tass(fs,mlop);
  typedef Dune::PDELab::OneStepGlobalAssembler<SASS,TASS,false> ASSEMBLER;
  ASSEMBLER assembler(sass,tass);

  // make a degree of freedom vector and set initial value
  typedef typename FS::DOF V;
  V x(fs.getGFS(),0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid.leafGridView(),problem);
  problem.setTime(0.0);
  Dune::PDELab::interpolate (g,fs.getGFS(),x);

  // linear solver backend
  typedef Dune::PDELab::ISTLSolverBackend_ExplicitDiagonal<FS,ASSEMBLER,solvertype> SBE;
  SBE sbe(fs,assembler,5000,1);

  // linear problem solver
  typedef Dune::PDELab::StationaryLinearProblemSolver<typename ASSEMBLER::GO,typename SBE::LS,V> PDESOLVER;
  PDESOLVER pdesolver(*assembler,*sbe,1e-6);

  // time-stepper
  //Dune::PDELab::OneStepThetaParameter<NumberType> method(0.5);
  //Dune::PDELab::Alexander3Parameter<NumberType> method;
  //typedef Dune::PDELab::OneStepMethod<NumberType,typename ASSEMBLER::GO,PDESOLVER,V> OSM;
  //OSM osm(method,*assembler,pdesolver);

  Dune::PDELab::ExplicitEulerParameter<NumberType> method1;
  Dune::PDELab::HeunParameter<NumberType> method2;
  Dune::PDELab::Shu3Parameter<NumberType> method3;
  Dune::PDELab::RK4Parameter<NumberType> method4;
  typedef Dune::PDELab::SimpleTimeController<NumberType> TC;
  TC tc;
  typedef Dune::PDELab::ExplicitOneStepMethod<NumberType,typename ASSEMBLER::GO,typename SBE::LS,V,V,TC> OSM;
  OSM osm(method2,*assembler,*sbe,tc);
  osm.setVerbosityLevel(2);

  // graphics for initial guess
  Dune::PDELab::FilenameHelper fn(basename);
  { // start a new block to automatically delete the VTKWriter object
    Dune::SubsamplingVTKWriter<typename GM::LeafGridView> vtkwriter(grid.leafGridView(),degree-1);
    typename FS::DGF xdgf(fs.getGFS(),x);
    vtkwriter.addVertexData(new typename FS::VTKF(xdgf,"x_h"));
    vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
    fn.increment();
  }

  // time loop
  NumberType time = 0.0;
  while (time<T-1e-10)
    {
      // assemble constraints for new time step (assumed to be constant for all substeps)
      problem.setTime(time+dt);
      fs.assembleConstraints(bctype);

      // do time step
      V xnew(fs.getGFS(),0.0);
      //osm.apply(time,dt,x,g,xnew);
      std::cout << "calling one step apply" << std::endl;
      osm.apply(time,dt,x,xnew);

      // output to VTK file
      {
        Dune::SubsamplingVTKWriter<typename GM::LeafGridView> vtkwriter(grid.leafGridView(),degree-1);
        typename FS::DGF xdgf(fs.getGFS(),xnew);
        vtkwriter.addVertexData(new typename FS::VTKF(xdgf,"x_h"));
        vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
        fn.increment();
      }

      // accept time step
      x = xnew;
      time += dt;
    }
}

//***********************************************************************
//***********************************************************************
// the main function
//***********************************************************************
//***********************************************************************

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // read command line arguments
  if (argc!=4)
    {
      std::cout << "usage: " << argv[0] << " <T> <dt> <cells>" << std::endl;
      return 0;
    }
  double T; sscanf(argv[1],"%lg",&T);
  double dt; sscanf(argv[2],"%lg",&dt);
  int cells; sscanf(argv[3],"%d",&cells);

  // start try/catch block to get error messages from dune
  try {

    const int dim=2;
    const int degree=1;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::sequential;
    const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
    const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;

    typedef Dune::YaspGrid<dim> GM;
    typedef Dune::PDELab::StructuredGrid<GM> Grid;
    Grid grid(elemtype,cells);
    grid->loadBalance();

    std::stringstream basename;
    basename << "advection_instationary" << "_dim" << dim << "_degree" << degree;
    do_simulation<GM,degree,elemtype,meshtype,solvertype>(T,dt,*grid,basename.str());
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
