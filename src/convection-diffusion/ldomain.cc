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
#include <dune/pdelab/gridfunctionspace/constraints.hh>
#include <dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/stationary/linearproblem.hh>


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

//! solve the diffusion problem on a given GridView with a given finite element space
template<class GV, class FEM, class CON>
void driver (const GV& gv, FEM& fem, std::string filename_base, int& N, double& l2e, double& h1se, double& ee)
{
  // coordinate and result type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // make grid function space 
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  N = gfs.globalSize();

  // make problem
  typedef ReentrantCornerProblem<GV,double> Problem;
  Problem problem;

  // make a std::vector of degree of freedom vectors and initialize it with Dirichlet extension
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0);
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
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // make linear solver and solve problem
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GFS> LS;
  LS ls (2, 5000, 1);
  // typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
  // LS ls(10000,1);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  SLP slp(gos,u,ls,1e-10);
  slp.apply();

  // compute errors
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  typedef DifferenceSquaredAdapter<G,DGF> DifferenceSquared;
  DifferenceSquared differencesquared(g,udgf);
  typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
  Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,10);
  l2e = sqrt(l2errorsquared);
  typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS,U> DGFGrad;
  DGFGrad udgfgrad(gfs,u);
  typedef ExactGradient<GV,double> Grad;
  Grad grad(gv);
  typedef DifferenceSquaredAdapter<Grad,DGFGrad> GradDifferenceSquared;
  GradDifferenceSquared graddifferencesquared(grad,udgfgrad);
  typename GradDifferenceSquared::Traits::RangeType h1semierrorsquared(0.0);
  Dune::PDELab::integrateGridFunction(graddifferencesquared,h1semierrorsquared,10);
  h1se = sqrt(h1semierrorsquared);

  // compute estimated error
  typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,Real,dim> P0FEM;
  P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::simplex,dim));
  typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,VBE> P0GFS; 
  P0GFS p0gfs(gv,p0fem);

  typedef Dune::PDELab::ConvectionDiffusionFEMResidualEstimator<Problem> ESTLOP;
  ESTLOP estlop(problem);
  typedef Dune::PDELab::GridOperatorSpace<GFS,P0GFS,ESTLOP,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,MBE> ESTGOS;
  ESTGOS estgos(gfs,p0gfs,estlop);

  typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,Real>::Type U0;
  U0 eta(p0gfs,0.0);
  estgos.residual(u,eta);
  ee = sqrt(eta.one_norm());

  // write vtk file
  // typedef DifferenceAdapter<G,DGF> Difference;
  // Difference difference(g,udgf);
  // Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
  // //Dune::VTKWriter<GV> vtkwriter(gv);
  // vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"u_h"));
  // vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<G>(g,"u"));
  // vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<Difference>(difference,"u-u_h"));
  // vtkwriter.write(filename_base,Dune::VTKOptions::binaryappended);
}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit 
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {

    // UGGrid/Pk2d version
#if HAVE_UG
    if (true)
    {
      // read command line arguments
      if (argc!=3)
	{
	  std::cout << "usage: " << argv[0] << " <mesh file> <maxlevel>" << std::endl;
	  return 0;
	}
      std::string grid_file(argv[1]);
      int max_level; sscanf(argv[2],"%d",&max_level);

      // make uggrid
      const int dim=2;
      typedef Dune::UGGrid<dim> GridType;
      GridType grid;
      typedef std::vector<int> GmshIndexMap;
      GmshIndexMap boundary_index_map;
      GmshIndexMap element_index_map;
      Dune::GmshReader<GridType> gmsh_reader;
      gmsh_reader.read(grid,grid_file,boundary_index_map,element_index_map,true,false);
      
      std::vector<double> l2(max_level+1,0.0);
      std::vector<double> h1s(max_level+1,0.0);
      std::vector<double> ee(max_level+1,0.0);
      std::vector<int> N(max_level+1,0);

      for (int i=0; i<=max_level; i++)
	{
	  typedef GridType::LeafGridView GV;
	  const GV& gv=grid.leafView();

	  // make finite element space and solve
	  typedef Dune::PDELab::ConformingDirichletConstraints CON;
	  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,GridType::ctype,double,1,dim> FEM;
	  FEM fem(gv);
	  std::stringstream filename;
	  filename << "ldomain_" << grid_file << "_l" << i;
	  driver<GV,FEM,CON>(gv,fem,filename.str(),N[i],l2[i],h1s[i],ee[i]);
	  if (i<max_level) grid.globalRefine(1);
	}

    std::cout << "Results on mesh=" << grid_file << std::endl; 
    std::cout << "N l2 l2rate h1semi h1semirate estimator effectivity" << std::endl;
    for (int i=0; i<=max_level; i++) 
      {
	double rate1=0.0;
	if (i>0) rate1=log(l2[i]/l2[i-1])/log(0.5);
	double rate2=0.0;
	if (i>0) rate2=log(h1s[i]/h1s[i-1])/log(0.5); 
	std::cout << std::setw(10) << N[i]
		  << std::setw(12) << std::setprecision(4) << std::scientific << l2[i]
		  << std::setw(12) << std::setprecision(4) << std::scientific << rate1
		  << std::setw(12) << std::setprecision(4) << std::scientific << h1s[i]
		  << std::setw(12) << std::setprecision(4) << std::scientific << rate2
		  << std::setw(12) << std::setprecision(4) << std::scientific << ee[i]
		  << std::setw(12) << std::setprecision(4) << std::scientific << ee[i]/(h1s[i])
		  << std::endl;
	}
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
