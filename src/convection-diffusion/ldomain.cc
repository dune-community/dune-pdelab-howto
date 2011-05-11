#include <iostream>
#include "config.h"           // file constructed by ./configure script
#include <dune/common/mpihelper.hh> // include mpi helper class 
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#if HAVE_UG 
#include <dune/grid/uggrid.hh>
#endif
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/constraints.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/adaptivity/adapt.hh>

#include"../utility/gridexamples.hh"

template<typename GV, typename RF>
class ReentrantCornerProblem
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
        I[i][j] = (i==j) ? 1.0 : 0.0;
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
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);

    typename Traits::DomainFieldType theta = std::atan2(x[1], x[0]);
    if(theta < 0.0) theta += 2*M_PI;
    typename Traits::DomainFieldType r = x.two_norm();

    return pow(r,2.0/3.0)*std::sin(theta*2.0/3.0);
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

//! exact gradient of solution
template<typename GV, typename RF>
class ExactGradient
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension,Dune::FieldVector<RF,GV::dimension> >,
      ExactGradient<GV,RF> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension,Dune::FieldVector<RF,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,ExactGradient<GV,RF> > BaseT;

  ExactGradient (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    typename Traits::DomainType x = e.geometry().global(xlocal);

    typename Traits::DomainFieldType theta = std::atan2(x[1], x[0]);
    if(theta < 0.0) theta += 2*M_PI;
    typename Traits::DomainFieldType r = x.two_norm();

    y[0] = (2.0/3.0)*pow(r,-1.0/3.0)*(cos(theta)*sin(2.0*theta/3.0) - sin(theta)*cos(2.0*theta/3.0));
    y[1] = (2.0/3.0)*pow(r,-1.0/3.0)*(sin(theta)*sin(2.0*theta/3.0) + cos(theta)*cos(2.0*theta/3.0));
  }
  
  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
};

/*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceAdapter 
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
				   typename T1::Traits::RangeFieldType,
				   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor 
  DifferenceAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x, 
                        typename Traits::RangeType& y) const
  {  
    typename Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }
  
private:
  const T1& t1;
  const T2& t2;
};

/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter 
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
				   typename T1::Traits::RangeFieldType,
				   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceSquaredAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor 
  DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x, 
                        typename Traits::RangeType& y) const
  {  
    typename T1::Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename T2::Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1.two_norm2();
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }
  
private:
  const T1& t1;
  const T2& t2;
};

//! Solve problem on leaf grid view and adapt grid
template<class Grid>
void driver (Grid& grid, std::string filename_base, double TOL, int maxsteps, double fraction)
{
  // some types
  typedef typename Grid::LeafGridView GV;
  typedef typename Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // some arrays to store results
  std::vector<double> l2;
  std::vector<double> h1s;
  std::vector<double> ee;
  std::vector<int> N;

  // make finite element map
  // note: adaptivity currently relies on finite element map not depending on grid view
  typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;

  // make grid function space 
  // note: adaptivity relies on leaf grid view object being updated by the grid on adaptation 
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(grid.leafView(),fem);
  N.push_back(gfs.globalSize());

  // make a degree of freedom vector;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0);

  // refinement loop
  for (int step=0; step<maxsteps; step++)
    {
      // get current leaf view
      const GV& gv=grid.leafView();

      // make problem
      typedef ReentrantCornerProblem<GV,Real> Problem;
      Problem problem;

      // initialize DOFS from Dirichlet extension
      // note: currently we start from scratch on every grid
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
      G g(gv,problem);
      Dune::PDELab::interpolate(g,gfs,u);

      // make constraints container and initialize it
      typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;
      Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv,problem);
      Dune::PDELab::constraints(bctype,gfs,cc);
      Dune::PDELab::set_nonconstrained_dofs(cc,0.0,u);
      std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;

      // make local operator
      typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
      LOP lop(problem);
      typedef VBE::MatrixBackend MBE;
      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
      GO go(gfs,cc,gfs,cc,lop);

      // make linear solver and solve problem
      typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
      LS ls (5000,1);
      // typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
      // LS ls(10000,1);
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      SLP slp(go,u,ls,1e-10);
      slp.apply();

      // compute errors
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
      DGF udgf(gfs,u);
      typedef DifferenceSquaredAdapter<G,DGF> DifferenceSquared;
      DifferenceSquared differencesquared(g,udgf);
      typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
      Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,10);
      l2.push_back(sqrt(l2errorsquared));
      typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,U> DGFGrad;
      DGFGrad udgfgrad(gfs,u);
      typedef ExactGradient<GV,Real> Grad;
      Grad grad(gv);
      typedef DifferenceSquaredAdapter<Grad,DGFGrad> GradDifferenceSquared;
      GradDifferenceSquared graddifferencesquared(grad,udgfgrad);
      typename GradDifferenceSquared::Traits::RangeType h1semierrorsquared(0.0);
      Dune::PDELab::integrateGridFunction(graddifferencesquared,h1semierrorsquared,10);
      h1s.push_back(sqrt(h1semierrorsquared));
      
      // compute estimated error
      typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> P0FEM;
      P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
      typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,VBE> P0GFS; 
      P0GFS p0gfs(gv,p0fem);
      typedef Dune::PDELab::ConvectionDiffusionFEMResidualEstimator<Problem> ESTLOP;
      ESTLOP estlop(problem);
      typedef Dune::PDELab::EmptyTransformation NoTrafo;
      typedef Dune::PDELab::GridOperator<GFS,P0GFS,ESTLOP,MBE,Real,Real,Real,NoTrafo,NoTrafo> ESTGO;
      ESTGO estgo(gfs,p0gfs,estlop);
      typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,Real>::Type U0;
      U0 eta(p0gfs,0.0);
      estgo.residual(u,eta);
      for (unsigned int i=0; i<eta.N(); i++) eta[i] = sqrt(eta[i]); 
      Real estimated_error = eta.two_norm();
      ee.push_back(estimated_error);

      // write vtk file
      typedef Dune::PDELab::DiscreteGridFunction<P0GFS,U0> DGF0;
      DGF0 udgf0(p0gfs,eta);
      typedef DifferenceAdapter<G,DGF> Difference;
      Difference difference(g,udgf);
      //Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
      Dune::VTKWriter<GV> vtkwriter(gv);
      std::stringstream fullname;
      fullname << filename_base << "_step" << step;
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"u_h"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<G>(g,"u"));
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<Difference>(difference,"u-u_h"));
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(udgf0,"estimated error"));
      vtkwriter.write(fullname.str(),Dune::VTKOptions::binaryappended);

      // error control
      if (estimated_error <= TOL) break;

      // create grid adaption classes
      typedef Dune::PDELab::L2Projection<GFS,U> Proj;
      Proj proj;
      typedef Dune::PDELab::ResidualErrorEstimation<GFS,U,ESTLOP> REE;
      REE ree(gfs,estlop);
      typedef Dune::PDELab::EstimationAdaptation<Grid,GFS,U,REE> EA;
      EA ea(grid,gfs,ree,fraction,0.0);
      typedef Dune::PDELab::GridAdaptor<Grid,GFS,U,EA,Proj> GRA;
      GRA gra(grid,gfs,ea,proj);
      
      // adapt grid
      if (step<maxsteps-1) 
	{
	  gra.adapt(u);
	  N.push_back(gfs.globalSize());
	}
    }

  // print results
  std::cout << "Results on mesh=" << filename_base << std::endl; 
  std::cout << "N l2 l2rate h1semi h1semirate estimator effectivity" << std::endl;
  for (std::size_t i=0; i<N.size(); i++) 
    {
      double rate1=0.0;
      if (i>0) rate1=log(l2[i]/l2[i-1])/log(0.5);
      double rate2=0.0;
      if (i>0) rate2=log(h1s[i]/h1s[i-1])/log(0.5); 
      std::cout << std::setw(10) << i
		<< std::setw(10) << N[i]
		<< std::setw(12) << std::setprecision(4) << std::scientific << l2[i]
		<< std::setw(12) << std::setprecision(4) << std::scientific << rate1
		<< std::setw(12) << std::setprecision(4) << std::scientific << h1s[i]
		<< std::setw(12) << std::setprecision(4) << std::scientific << rate2
		<< std::setw(12) << std::setprecision(4) << std::scientific << ee[i]
		<< std::setw(12) << std::setprecision(4) << std::scientific << ee[i]/(h1s[i])
		<< std::endl;
    }
  

}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit 
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {

    // Alberta/Pk2d version
#if HAVE_ALBERTA
    if (true)
    {
      // read command line arguments
      if (argc!=6)
	{
	  std::cout << "usage: " << argv[0] << " <mesh file> <start level> <TOL> <max steps> <fraction>" << std::endl;
	  return 0;
	}
      std::string grid_file(argv[1]);
      int start_level; sscanf(argv[2],"%d",&start_level);
      double TOL; sscanf(argv[3],"%lg",&TOL);
      int maxsteps; sscanf(argv[4],"%d",&maxsteps);
      double fraction; sscanf(argv[5],"%lg",&fraction);

      // make Alberta grid
      typedef AlbertaLDomain::Grid GridType;
      AlbertaLDomain gridp;
      GridType &grid = gridp;

      // make UG grid
      // const int dim=2;
      // typedef Dune::UGGrid<dim> GridType;
      // GridType grid;
      // typedef std::vector<int> GmshIndexMap;
      // GmshIndexMap boundary_index_map;
      // GmshIndexMap element_index_map;
      // Dune::GmshReader<GridType> gmsh_reader;
      // gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);

      for (int i=0; i<start_level; i++) grid.globalRefine(1);
      driver(grid,"ldomain",TOL,maxsteps,fraction);
    }
#endif
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
