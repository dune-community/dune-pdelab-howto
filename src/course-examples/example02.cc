// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem in constrained spaces with conforming finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<map>
#include<string>
#include<vector>

#include<math.h>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>

#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif

#include<dune/istl/bvector.hh>
#include<dune/istl/io.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"example02_bctype.hh"
#include"example02_bcextension.hh"
#include"example02_operator.hh"
#include"example02_Q1.hh"

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

    if (argc!=2)
      {
        if(helper.rank()==0)
          std::cout << "usage: ./example02 <level>" << std::endl;
        return 1;
      }

    int level;
    sscanf(argv[1],"%d",&level);

    // sequential version
    if (1 && helper.size()==1)
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      example02_Q1(gv);
    }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    throw;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    throw;
  }
}
