// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    
    \brief Solve parabolic problem with conforming finite elements in space and implicit one step methods in time
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG 
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include"example03_bctype.hh"
#include"example03_bcextension.hh"
#include"example03_operator.hh"
#include"example03_toperator.hh"
#include"example03_Q2.hh"

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

    if (argc!=4)
      {
        if(helper.rank()==0) {
          std::cout << "usage:      ./example03 <level> <dt> <tend>" << std::endl;
          std::cout << "suggestion: level=3 dt=0.02 tend=6.28" << std::endl;
        }
        return 1;
      }

    int level;
    sscanf(argv[1],"%d",&level);

    double dt;
    sscanf(argv[2],"%lg",&dt);

    double tend;
    sscanf(argv[3],"%lg",&tend);

    // sequential version
    if (helper.size()==1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      example03_Q2(gv,dt,tend);
    }
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
