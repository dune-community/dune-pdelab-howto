// -*- tab-width: 4; indent-tabs-mode: nil -*-

#include <iostream>
#include <sstream>
#include "config.h"           // file constructed by ./configure script
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
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
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/finiteelementmap/opbfem.hh>
#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include<dune/pdelab/localoperator/errorindicatordg.hh>

#include<dune/pdelab/stationary/linearproblem.hh>

#include<dune/pdelab/adaptivity/adaptivity.hh>

#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>

#include"../utility/gridexamples.hh"

#include "reentrantcornerproblem.hh"

//! Solve problem on leaf grid view and adapt grid
template<typename Grid,typename FEMDG,int degree,int blocksize>
void driverDG ( Grid& grid,
                const Dune::GeometryType& gt,
                const FEMDG& femdg,
                const Dune::ParameterTree& configuration )
{

  int verbose = configuration.get<int>("general.verbose");

  int strategy = configuration.get<int>("adaptivity.strategy");
  std::stringstream vtu,cmd;
  vtu << "dgldomain_s" << strategy;
  std::string filename_base (vtu.str());
  cmd << "rm dgldomain_s" << strategy << "*.vtu";
  system( cmd.str().c_str() );

  // some types
  typedef typename Grid::LeafGridView GV;
  typedef typename Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  typedef ReentrantCornerProblem<GV,Real> Problem;
  Problem problem;
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> ESOL;

  // some arrays to store results
  std::vector<double> l2;
  std::vector<double> h1s;
  std::vector<double> ee;
  std::vector<int> N;
  std::vector<int> nIterations; // of the linear solver

  // make grid function space
  // note: adaptivity relies on leaf grid view object being updated by the grid on adaptation
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> VBE1;
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,blocksize> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE> GFS;
  GFS gfs(grid.leafGridView(),femdg);
  gfs.update();
  N.push_back(gfs.globalSize());

  // some local operator parameters
  double alpha = 2.0;
  Dune::PDELab::ConvectionDiffusionDGMethod::Type m
    = Dune::PDELab::ConvectionDiffusionDGMethod::SIPG;
  Dune::PDELab::ConvectionDiffusionDGWeights::Type w
    = Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn;


  // make a degree of freedom vector;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;


  int maxsteps = configuration.get<int>("adaptivity.maxsteps");
  bool bReplay = configuration.get<int>("adaptivity.replay");
  double TOL = configuration.get<double>("adaptivity.TOL");
  double refinementfraction = configuration.get<double>("adaptivity.refinementfraction");

  // refinement loop
  for (int step=0; step<maxsteps; step++)
    {

      std::cout << "***************************************" << std::endl;
      std::cout << "Refinement Step " << step << std::endl;
      std::cout << "***************************************" << std::endl;

      // get current leaf view
      const GV& gv=grid.leafGridView();

      // make constraints container and initialize it
      typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
      CC cc;

      if( bReplay==false ) {
        // reset degree of freedom vector;
        std::cout << "Resetting solution on the refined grid." << std::endl;
        u = U(gfs,0.0);
      }
      else {
        std::cout << "Solution from coarse grid transferred to refined grid." << std::endl;
        // write vtk file
        DGF pre_udgf(gfs,u);
        Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,degree-1);
        std::stringstream fullname;
        fullname << filename_base << "_prestep" << step;
        // plot analytical solution on the refined gridview
        vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(pre_udgf,"pre_u_h"));
        vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);
      }

      // make local operator
      typedef Dune::PDELab::ConvectionDiffusionDG<Problem,FEMDG> LOP;
      LOP lop(problem,m,w,alpha);
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
      GO go(gfs,cc,gfs,cc,lop,mbe);

      // make linear solver and solve problem
      //typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
      //LS ls (gfs,50,2);
      typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
      LS ls(10000,verbose);
      typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
      double reduction = configuration.get<double>("istl.reduction");
      SLP slp(go,ls,u,reduction);
      slp.apply();
      nIterations.push_back( slp.ls_result().iterations );


      // compute errors
      ESOL exactsolution(gv,problem);
      DGF udgf(gfs,u);
      typedef DifferenceSquaredAdapter<ESOL,DGF> DifferenceSquared;
      DifferenceSquared differencesquared(exactsolution,udgf);
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

      P0FEM p0fem(gt);

      typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,Dune::PDELab::NoConstraints,VBE1> P0GFS;
      P0GFS p0gfs(gv,p0fem);

      //typedef Dune::PDELab::ConvectionDiffusionFEMResidualEstimator<Problem> ESTLOP;
      //ESTLOP estlop(problem);

      typedef Dune::PDELab::ConvectionDiffusionDG_ErrorIndicator<Problem> ESTLOP;
      ESTLOP estlop(problem,m,w,alpha);

      typedef Dune::PDELab::EmptyTransformation NoTrafo;
      typedef Dune::PDELab::GridOperator<GFS,P0GFS,ESTLOP,MBE,Real,Real,Real,NoTrafo,NoTrafo> ESTGO;
      ESTGO estgo(gfs,p0gfs,estlop,mbe);
      typedef typename Dune::PDELab::BackendVectorSelector<P0GFS,Real>::Type U0;
      U0 eta(p0gfs,0.0);
      estgo.residual(u,eta);

      for( typename U0::iterator it = eta.begin(), end = eta.end();
           it != end; ++it )
        *it = sqrt(*it);

      Real estimated_error = eta.two_norm();
      ee.push_back(estimated_error);

      // write vtk file
      typedef Dune::PDELab::DiscreteGridFunction<P0GFS,U0> DGF0;
      DGF0 udgf0(p0gfs,eta);
      typedef DifferenceAdapter<ESOL,DGF> Difference;
      Difference difference(exactsolution,udgf);
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,degree-1);
      //Dune::VTKWriter<GV> vtkwriter(gv);
      std::stringstream fullname;
      fullname << filename_base << "_step" << step;
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"u_h"));
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<ESOL>(exactsolution,"u"));
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<Difference>(difference,"u-u_h"));
      vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF0>(udgf0,"estimated error"));
      vtkwriter.write(fullname.str(),Dune::VTK::appendedraw);

      // error control
      if (estimated_error <= TOL) break;

      // adapt grid
      if (step<maxsteps-1) {
        double alpha_r(refinementfraction);
        double refine_threshold(0);
        double beta_c(0);       // no coarsening here
        double eta_beta(0);   // dummy
        if( strategy == 1)
          error_fraction( eta, alpha_r, beta_c, refine_threshold, eta_beta, verbose );
        else
          element_fraction( eta, alpha_r, beta_c, refine_threshold, eta_beta, verbose );
        mark_grid( grid, eta, refine_threshold, 0.0 );
        adapt_grid( grid, gfs, u, 2 * degree );
        N.push_back( gfs.globalSize() );
      }

    }

  // print results
  std::cout << "Results for DG polynomial degree " << degree
            << " with blocksize " << blocksize
            << std::endl;
  std::cout << "           N"
            << "    IT"
            << "          l2"
            << "      l2rate"
            << "      h1semi"
            << "  h1semirate"
            << "   estimator"
            << " effectivity" << std::endl;
  for (std::size_t i=0; i<N.size(); i++)
    {
      double rate1=0.0;
      if (i>0) rate1=log(l2[i]/l2[i-1])/log(0.5);
      double rate2=0.0;
      if (i>0) rate2=log(h1s[i]/h1s[i-1])/log(0.5);
      std::cout << std::setw(3) << i
                << std::setw(9) << N[i]
                << std::setw(6) << nIterations[i]
                << std::setw(12) << std::setprecision(4) << std::scientific << l2[i]
                << std::setw(12) << std::setprecision(4) << std::scientific << rate1
                << std::setw(12) << std::setprecision(4) << std::scientific << h1s[i]
                << std::setw(12) << std::setprecision(4) << std::scientific << rate2
                << std::setw(12) << std::setprecision(4) << std::scientific << ee[i]
                << std::setw(12) << std::setprecision(4) << std::scientific << ee[i]/(h1s[i])
                << std::endl;
    }

  std::cout << "View results using: \n paraview --data=" << vtu.str() << "_step..vtu" << std::endl;

}

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try {

    std::string config_file("ldomain.ini");
    Dune::ParameterTree configuration;
    Dune::ParameterTreeParser parser;
    try{
      parser.readINITree( config_file, configuration );
    }
    catch(...){
      std::cerr << "Could not read config file \""
                << config_file << "\"!" << std::endl;
      exit(1);
    }

    // make finite element map
    // note: adaptivity currently relies on finite element map
    // not depending on grid view

    // Hint: Set the polynomial degree here
    const int degree=1;

    // UG version
#if HAVE_UG
    if( "ug"==configuration.get<std::string>("grid.manager") ) {
      // make UG grid
      const int dim=2;

      if( "cube"==configuration.get<std::string>("ug.geometrytype") ){
        // example of a simple mesh based on cubes
        typedef UGLDomainCubes GridType;
        GridType grid(1000);
        std::cout << "cubes: only non-conforming refinement possible here" << std::endl;
        grid.setClosureType( Dune::UGGrid<dim>::NONE );
        grid.globalRefine( configuration.get<int>("grid.baselevel") );
        //grid.loadBalance();
        Dune::GeometryType gt;
        gt.makeCube(dim);
        //typedef Dune::PDELab::OPBLocalFiniteElementMap<GridType::ctype,double,degree,dim,Dune::GeometryType::cube> FEMDG;
        typedef Dune::PDELab::QkDGLocalFiniteElementMap<GridType::ctype,double,degree,dim> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::QkStuff::QkSize<degree,dim>::value;
        driverDG<GridType,FEMDG,degree,blocksize>(grid,gt,femdg,configuration);
      }
      else{
        /*
          // example of a mesh generated by gmsh:
          typedef Dune::UGGrid<dim> GridType;
          GridType grid;
          typedef std::vector<int> GmshIndexMap;
          GmshIndexMap boundary_index_map;
          GmshIndexMap element_index_map;
          Dune::GridFactory<GridType> factory(&grid);
          std::string grid_file="grids/ldomain.msh";
          Dune::GmshReader<GridType>::read(factory,grid_file,boundary_index_map,element_index_map,true,false);
          factory.createGrid();
        */

        // example of a simple mesh based on simplices
        typedef UGLDomain GridType;
        GridType grid(1000);
        if( "nonconforming"==configuration.get<std::string>("ug.refinementtype") ){
          std::cout << "non-conforming refinement" << std::endl;
          grid.setClosureType( Dune::UGGrid<dim>::NONE );
        }
        grid.globalRefine( configuration.get<int>("grid.baselevel") );
        //grid.loadBalance();
        Dune::GeometryType gt;
        gt.makeSimplex(dim);
        typedef Dune::PDELab::OPBLocalFiniteElementMap<GridType::ctype,double,degree,dim,Dune::GeometryType::simplex> FEMDG;
        FEMDG femdg;
        const int blocksize = Dune::PB::PkSize<degree,dim>::value;
        driverDG<GridType,FEMDG,degree,blocksize>(grid,gt,femdg,configuration);
      }

    }
#endif

#if HAVE_ALBERTA
    if( "alberta"==configuration.get<std::string>("grid.manager") ) {
      // make Alberta grid
      const int dim=2;
      typedef AlbertaLDomain::Grid GridType;
      AlbertaLDomain gridp;
      GridType &grid = gridp;
      grid.globalRefine( configuration.get<int>("grid.baselevel") );

      // Finite Element Map
      typedef Dune::PDELab::MonomLocalFiniteElementMap<GridType::ctype,double,dim,degree> FEMDG;
      Dune::GeometryType gt;
      gt.makeSimplex(dim);
      FEMDG femdg(gt);
      const int blocksize = Dune::PB::PkSize<degree,dim>::value;
      driverDG<GridType,FEMDG,degree,blocksize>(grid,gt,femdg,configuration);
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
