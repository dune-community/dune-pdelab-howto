// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

/** Parameter class for the stationary convection-diffusion equation of the following form:
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \ \
 *                                              u &=& g \mbox{ on } \partial\Omega_D (Dirichlet)\ \
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N (Flux)\ \
 *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O (Outflow)
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * The template parameters are:
 *  - GV a model of a GridView
 *  - RF numeric type to represent results
 */

template<typename GV, typename RF>
class GenericEllipticProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
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
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);
    return exp(-(x*x));
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
};

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // command line args
  int cells=10; if (argc>=2) sscanf(argv[1],"%d",&cells);

  // define parameters
  const unsigned int dim = 3;
  const unsigned int degree = 2;
  const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
  const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::sequential;
  typedef double NumberType;

  // make grid
#if HAVE_UG
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
  typedef Dune::UGGrid<dim> GM;
#elif HAVE_ALUGRID
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
  typedef Dune::ALUGrid<dim,dim,Dune::simplex,Dune::nonconforming> GM;
#else  // ! (HAVE_UG || HAVE_ALUGRID)
  const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
  typedef Dune::YaspGrid<dim> GM;
#endif
  typedef Dune::PDELab::StructuredGrid<GM> Grid;
  Grid grid(elemtype,cells);
  grid->loadBalance();

  // make problem parameters
  typedef GenericEllipticProblem<GM::LeafGridView,NumberType> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
  BCType bctype(grid->leafGridView(),problem);

  // make a finite element space
  typedef Dune::PDELab::CGSpace<GM,NumberType,degree,BCType,elemtype,meshtype,solvertype> FS;
  FS fs(*grid,bctype);

  // make a degree of freedom vector and initialize it with a function
  typedef FS::DOF V;
  V x(fs.getGFS(),0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
  G g(grid->leafGridView(),problem);
  Dune::PDELab::interpolate(g,fs.getGFS(),x);

  // assemble constraints
  fs.assembleConstraints(bctype);
  fs.setNonConstrainedDOFS(x,0.0);

  // assembler for finite elemenent problem
  typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FS::FEM> LOP;
  LOP lop(problem);
  typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASSEMBLER;
  ASSEMBLER assembler(fs,lop);

  // make linear solver and solve problem
  //typedef Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASSEMBLER,solvertype> SBE;
  typedef Dune::PDELab::ISTLSolverBackend_CG_AMG_SSOR<FS,ASSEMBLER,solvertype> SBE;
  SBE sbe(fs,assembler,5000,1);
  typedef Dune::PDELab::StationaryLinearProblemSolver<ASSEMBLER::GO,SBE::LS,V> SLP;
  SLP slp(*assembler,*sbe,x,1e-6);
  slp.apply();

  // output grid to VTK file
  Dune::SubsamplingVTKWriter<GM::LeafGridView> vtkwriter(grid->leafGridView(),degree-1);
  FS::DGF xdgf(fs.getGFS(),x);
  vtkwriter.addVertexData(new FS::VTKF(xdgf,"x_h"));
  vtkwriter.write("poisson_uniform",Dune::VTK::appendedraw);

  // done
  return 0;
}
