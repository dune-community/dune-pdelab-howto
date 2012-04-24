// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Solve Problems A-F in parallel on overlapping grids using conforming linear finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<string>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>

#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"../utility/gridexamples.hh"
#define PROBLEM_A

#ifdef PROBLEM_A
#include "parameterA.hh"
#endif

#ifdef PROBLEM_B
#include "parameterB.hh"
#endif

#ifdef PROBLEM_C
#include "parameterC.hh"
#endif

#ifdef PROBLEM_D
#include "parameterD.hh"
#endif

#ifdef PROBLEM_E
#include "parameterE.hh"
#endif

#ifdef PROBLEM_F
#include "parameterF.hh"
#endif

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template<typename PROBLEM, typename GV, typename FEM> 
void driver (PROBLEM& problem,
             const GV& gv, const FEM& fem, std::string filename)
{
  // constants and types and global variables
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type CC;
  CC cc;
  cc.clear();
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<PROBLEM> bctype(gv,problem);

  // make grid operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<PROBLEM,FEM> LOP; 
  LOP lop(problem);
  typedef Dune::PDELab::GridOperator<GFS,GFS,
    LOP,VBE::MatrixBackend,R,R,R,CC,CC> GO;           
  GO go(gfs,cc,gfs,cc,lop);                     

  // make coefficent Vector and initialize it from a function
  typedef typename GO::Traits::Domain V;
  V x(gfs,0.0);
  typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<PROBLEM> G;
  G g(gv,problem);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::constraints(bctype,gfs,cc);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,x);

  typedef  Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GO> LS;
  LS ls(gfs,5000,3);
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
  SLP slp(go,x,ls,1e-12);
  slp.apply();

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.pwrite(filename.c_str(),"vtk","",Dune::VTKOptions::binaryappended);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      {
        if(helper.rank()==0)
          std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
      }
    
    // Q1, 2d
    if (true)
      {
        // make grid
        Dune::FieldVector<double,2> L(1.0);
        Dune::FieldVector<int,2> N(3);
        Dune::FieldVector<bool,2> B(false);
        int overlap=1;
        Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
        //grid.globalRefine(4);
        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif
      
        typedef Dune::YaspGrid<2>::ctype DF;
        typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
        FEM fem;

        Dune::FieldVector<double,2> correlation_length;
        correlation_length = 1.0/64.0;
        driver( problem,
                gv,
                fem,
                "single_phase_yasp2d_Q1" );
      }

    // Q1, 3d
    if (false)
      {
        // make grid
        Dune::FieldVector<double,3> L(1.0);
        Dune::FieldVector<int,3> N(4);
        Dune::FieldVector<bool,3> B(false);
        int overlap=0;
        Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
        grid.globalRefine(4);
        typedef Dune::YaspGrid<3>::LeafGridView GV;
        const GV& gv=grid.leafView();
#ifdef PROBLEM_A
        typedef ParameterA<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_B
        typedef ParameterB<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_C
        typedef ParameterC<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_D
        typedef ParameterD<GV,double> PROBLEM;
        Dune::FieldVector<double,GV::Grid::dimension> correlation_length;
        correlation_length = 1.0/64.0;
        PROBLEM problem(gv,correlation_length,1.0,0.0,5000,-1083);
#endif
#ifdef PROBLEM_E
        typedef ParameterE<GV,double> PROBLEM;
        PROBLEM problem;
#endif
#ifdef PROBLEM_F
        typedef ParameterF<GV,double> PROBLEM;
        PROBLEM problem;
#endif

        typedef Dune::YaspGrid<3>::ctype DF;
        typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
        FEM fem;

        Dune::FieldVector<double,3> correlation_length;
        correlation_length = 1.0/64.0;
        driver( problem,
                gv,
                fem,
                "single_phase_yasp3d_Q1" );
      }

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
