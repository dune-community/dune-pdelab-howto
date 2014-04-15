// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>

//***********************************************************************
//***********************************************************************
// define the reentrant corner in the L-domain with known exact solution
//***********************************************************************
//***********************************************************************

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

//***********************************************************************
//***********************************************************************
// a grid function giving the gradient of the exact solution
// needed for computing H^1 errors
//***********************************************************************
//***********************************************************************

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

//***********************************************************************
//***********************************************************************
// some function adapters
//***********************************************************************
//***********************************************************************

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

//***********************************************************************
//***********************************************************************
// meshes for the ldomain
//***********************************************************************
//***********************************************************************

template<class Factory>
void ldomain_mesh (Factory& factory, Dune::GeometryType::BasicType elemtype)
{
  if (elemtype == Dune::GeometryType::simplex)
    {
      /* 5 - 6 - 7
         |   |   |
         2 - 3 - 4
         |   |
         0 - 1
      */
      Dune::FieldVector< double, 2 > pos;
      pos[0] = -1;  pos[1] = -1; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] = -1; factory.insertVertex(pos);
      pos[0] = -1;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] =  1;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] = -1;  pos[1] =  1; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] =  1; factory.insertVertex(pos);
      pos[0] =  1;  pos[1] =  1; factory.insertVertex(pos);

      const Dune::GeometryType type( Dune::GeometryType::simplex, 2 );
      std::vector< unsigned int > cornerIDs( 3 );
      cornerIDs[0] = 3; cornerIDs[1] = 0; cornerIDs[2] = 1;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 0; cornerIDs[1] = 3; cornerIDs[2] = 2;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 3; cornerIDs[1] = 5; cornerIDs[2] = 2;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 5; cornerIDs[1] = 3; cornerIDs[2] = 6;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 7; cornerIDs[1] = 3; cornerIDs[2] = 4;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 3; cornerIDs[1] = 7; cornerIDs[2] = 6;
      factory.insertElement( type, cornerIDs );
    }
  if (elemtype == Dune::GeometryType::cube)
    {
      /* 5 - 6 - 7
         |   |   |
         2 - 3 - 4
         |   |
         0 - 1
      */
      Dune::FieldVector< double, 2 > pos;
      pos[0] = -1;  pos[1] = -1; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] = -1; factory.insertVertex(pos);
      pos[0] = -1;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] =  1;  pos[1] =  0; factory.insertVertex(pos);
      pos[0] = -1;  pos[1] =  1; factory.insertVertex(pos);
      pos[0] =  0;  pos[1] =  1; factory.insertVertex(pos);
      pos[0] =  1;  pos[1] =  1; factory.insertVertex(pos);

      const Dune::GeometryType type( Dune::GeometryType::cube, 2 );
      std::vector< unsigned int > cornerIDs( 4 );
      cornerIDs[0] = 0; cornerIDs[1] = 1; cornerIDs[2] = 2; cornerIDs[3] = 3;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 2; cornerIDs[1] = 3; cornerIDs[2] = 5; cornerIDs[3] = 6;
      factory.insertElement( type, cornerIDs );
      cornerIDs[0] = 3; cornerIDs[1] = 4; cornerIDs[2] = 6; cornerIDs[3] = 7;
      factory.insertElement( type, cornerIDs );
    }
}

template<typename GM, unsigned int degree, Dune::GeometryType::BasicType elemtype,
         Dune::PDELab::MeshType meshtype, Dune::SolverCategory::Category solvertype>
void driver (Dune::shared_ptr<GM> grid, int prerefine_level, double tol, int maxsteps,
             double fraction, std::string basename)
{
  // define parameters
  typedef double NumberType;

  // make grid
  for (int i=0; i<prerefine_level; i++) grid->globalRefine(1);

  // some arrays to store result of adaptive computation
  std::vector<double> l2;
  std::vector<double> h1s;
  std::vector<double> ee;
  std::vector<int> N;
  std::vector<int> maxlevel;

  // adaptive loop
  for (int step=1; step<=maxsteps; step++)
    {
      std::cout << "*** adaptive step #" << step << std::endl;

      // make problem parameters
      typedef ReentrantCornerProblem<typename GM::LeafGridView,NumberType> Problem;
      Problem problem;
      typedef Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> BCType;
      BCType bctype(grid->leafGridView(),problem);

      // make a finite element space
      typedef Dune::PDELab::CGSpace<GM,NumberType,degree,BCType,elemtype,meshtype,solvertype> FS;
      FS fs(*grid,bctype);
      N.push_back(fs.getGFS().globalSize());
      maxlevel.push_back(grid->maxLevel());

      // make a degree of freedom vector
      typedef typename FS::DOF X;
      X x(fs.getGFS(),0.0);

      // assemble constraints and Dirichlet BC
      fs.assembleConstraints(bctype);
      typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
      G g(grid->leafGridView(),problem);
      Dune::PDELab::interpolate(g,fs.getGFS(),x);

      // assembler for finite elemenent problem
      typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,typename FS::FEM> LOP;
      LOP lop(problem);
      typedef Dune::PDELab::GalerkinGlobalAssembler<FS,LOP,solvertype> ASS;
      ASS ass(fs,lop);

      // make linear solver and solve problem
      // typedef Dune::PDELab::ISTLSolverBackend_CG_AMG_SSOR<FS,ASS,solvertype> SBEa;
      typedef Dune::PDELab::ISTLSolverBackend_IterativeDefault<FS,ASS,solvertype> SBE;
      SBE sbe(fs,ass,5000,1);

      // solve linear system in case of hanging nodes (this should go to StationaryLinearProblemSolver)
      typedef Dune::PDELab::StationaryLinearProblemSolver<typename ASS::GO,typename SBE::LS,X> SLP;
      SLP slp(*ass,*sbe,x,1e-8);
      slp.apply();

      // compute errors
      typename FS::DGF xdgf(fs.getGFS(),x);
      typedef DifferenceSquaredAdapter<G,typename FS::DGF> DifferenceSquared;
      DifferenceSquared differencesquared(g,xdgf);
      typename DifferenceSquared::Traits::RangeType l2errorsquared(0.0);
      Dune::PDELab::integrateGridFunction(differencesquared,l2errorsquared,8);
      l2.push_back(sqrt(l2errorsquared));
      typedef Dune::PDELab::DiscreteGridFunctionGradient<typename FS::GFS,X> DGFGrad;
      DGFGrad xdgfgrad(fs.getGFS(),x);
      typedef ExactGradient<typename GM::LeafGridView,NumberType> Grad;
      Grad grad(grid->leafGridView());
      typedef DifferenceSquaredAdapter<Grad,DGFGrad> GradDifferenceSquared;
      GradDifferenceSquared graddifferencesquared(grad,xdgfgrad);
      typename GradDifferenceSquared::Traits::RangeType h1semierrorsquared(0.0);
      Dune::PDELab::integrateGridFunction(graddifferencesquared,h1semierrorsquared,8);
      h1s.push_back(sqrt(h1semierrorsquared));

      // a posteriori error estimate
      typedef Dune::PDELab::P0Space<GM,NumberType,elemtype,solvertype> ESTFS;
      ESTFS estfs(grid->leafGridView());
      typedef Dune::PDELab::ConvectionDiffusionFEMResidualEstimator<Problem> ESTLOP;
      ESTLOP estlop(problem);
      typedef Dune::PDELab::GlobalAssembler<FS,ESTFS,ESTLOP,solvertype> ESTASS;
      ESTASS estass(fs,estfs,estlop);
      typedef typename ESTFS::DOF Y;
      Y y(estfs.getGFS(),0.0);
      estass->residual(x,y);
      NumberType estimated_error = sqrt(y.one_norm());
      std::cout << "estimated error = " << estimated_error << std::endl;
      ee.push_back(estimated_error);

      // output grid to VTK file
      Dune::VTKWriter<typename GM::LeafGridView> vtkwriter(grid->leafGridView());
      vtkwriter.addVertexData(new typename FS::VTKF(xdgf,"x_h"));
      Y ee(estfs.getGFS(),0.0);
      for (unsigned int i=0; i<ee.base().N(); i++) (ee.base())[i] = sqrt((y.base())[i]);
      typename ESTFS::DGF eedgf(estfs.getGFS(),ee);
      vtkwriter.addCellData(new typename ESTFS::VTKF(eedgf,"estimated error"));
      std::stringstream fullname;
      fullname << basename << step;
      vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);

      // error control
      if (estimated_error<tol) break;

      // mark elements for refinement
      NumberType eta_refine, eta_coarsen;
      Dune::PDELab::error_fraction(y,fraction,0.0,eta_refine,eta_coarsen);
      Dune::PDELab::mark_grid(*grid,y,eta_refine,eta_coarsen);

      // do refinement
      Dune::PDELab::adapt_grid(*grid,fs.getGFS(),x,2*degree);
    }

  // print results
  //  std::cout << "Results on mesh=" << filename_base << std::endl;
  std::cout << "#step N maxlevel l2 l2rate h1semi h1semirate estimator effectivity" << std::endl;
  for (std::size_t i=0; i<N.size(); i++)
    {
      double rate1=0.0;
      if (i>0) rate1=-log(l2[i]/l2[i-1])/log((1.0*N[i])/N[i-1]);
      double rate2=0.0;
      if (i>0) rate2=-log(h1s[i]/h1s[i-1])/log((1.0*N[i])/N[i-1]);
      std::cout << std::setw(10) << i
                << std::setw(10) << N[i]
                << std::setw(10) << maxlevel[i]
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

    // read command line arguments
    if (argc<7 || argc>8)
      {
        std::cout << "usage: " << argv[0] << " <mode> <degree> <start level> <TOL> <max steps> <fraction> [<basename>]" << std::endl;
        std::cout << "mode=1 Alberta" << std::endl;
        std::cout << "mode=2 UG conforming" << std::endl;
        std::cout << "mode=3 ALUSimplex" << std::endl;
        std::cout << "mode=4 UG quads nonconforming" << std::endl;
        std::cout << "degree polynomial degree for conforming calculation: 1,2,3,4" << std::endl;
        std::cout << "basename for output (optional)" << std::endl;
        return 0;
      }
    int mode; sscanf(argv[1],"%d",&mode);
    int degree; sscanf(argv[2],"%d",&degree);
    int start_level; sscanf(argv[3],"%d",&start_level);
    double TOL; sscanf(argv[4],"%lg",&TOL);
    int maxsteps; sscanf(argv[5],"%d",&maxsteps);
    double fraction; sscanf(argv[6],"%lg",&fraction);
    char bn[129];
    if (argc==8)
      sscanf(argv[7],"%s",bn);
    else
      sscanf("ldomain","%s",bn);

    const int dim=2;
    const Dune::SolverCategory::Category solvertype = Dune::SolverCategory::sequential;

#if HAVE_ALBERTA
    if (mode==1)
      {
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
        const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
        typedef Dune::AlbertaGrid<dim,dim> GM;
        Dune::GridFactory<GM> factory;
        ldomain_mesh(factory,elemtype);
        Dune::shared_ptr<GM> grid(factory.createGrid());
        std::stringstream basename;
        basename << bn << "_alberta_P" << degree;

        if (degree==1) driver<GM,1,elemtype,meshtype,solvertype>(grid,dim*start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==2) driver<GM,2,elemtype,meshtype,solvertype>(grid,dim*start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==3) driver<GM,3,elemtype,meshtype,solvertype>(grid,dim*start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==4) driver<GM,4,elemtype,meshtype,solvertype>(grid,dim*start_level,TOL,maxsteps,fraction,basename.str());
      }
#endif

#if HAVE_UG
    if (mode==2)
      {
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
        const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::conforming;
        typedef Dune::UGGrid<dim> GM;
        GM::setDefaultHeapSize(1500);
        Dune::GridFactory<GM> factory;
        ldomain_mesh(factory,elemtype);
        Dune::shared_ptr<GM> grid(factory.createGrid());
        std::stringstream basename;
        basename << bn << "_ug_conforming_P" << degree;

        if (degree==1) driver<GM,1,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==2) driver<GM,2,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==3) driver<GM,3,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
        if (degree==4) driver<GM,4,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
      }
    if (mode==4)
      {
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::cube;
        const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::nonconforming;
        typedef Dune::UGGrid<dim> GM;
        GM::setDefaultHeapSize(1500);
        Dune::GridFactory<GM> factory;
        ldomain_mesh(factory,elemtype);
        Dune::shared_ptr<GM> grid(factory.createGrid());
        grid->setRefinementType( Dune::UGGrid<dim>::LOCAL );
        grid->setClosureType( Dune::UGGrid<dim>::NONE );
        std::stringstream basename;
        basename << bn << "_ug_conforming_Q1";

        if (degree==1) driver<GM,1,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
      }
#endif

#if HAVE_ALUGRID
    if (mode==3)
      {
        const Dune::GeometryType::BasicType elemtype = Dune::GeometryType::simplex;
        const Dune::PDELab::MeshType meshtype = Dune::PDELab::MeshType::nonconforming;
        typedef Dune::ALUGrid<dim,dim,Dune::simplex,Dune::nonconforming> GM;
        Dune::GridFactory<GM> factory;
        ldomain_mesh(factory,elemtype);
        Dune::shared_ptr<GM> grid(factory.createGrid());
        std::stringstream basename;
        basename << bn << "_alu_nonconforming_P1";

        if (degree==1) driver<GM,1,elemtype,meshtype,solvertype>(grid,start_level,TOL,maxsteps,fraction,basename.str());
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
