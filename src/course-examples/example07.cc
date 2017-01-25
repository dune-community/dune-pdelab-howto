// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem with adaptive conforming finite element method
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
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include<dune/pdelab/adaptivity/adaptivity.hh>

#include"example02_bctype.hh"
#include"example02_bcextension.hh"
#include"example02_operator.hh"
#include"example07_error_indicator.hh"
#include"example07_adaptivity.hh"

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

    if (argc!=3)
      {
        if(helper.rank()==0)
          std::cout << "usage: ./example07 <startLevel> <maxLevel>" << std::endl;
        return 1;
      }

    int startLevel;
    sscanf(argv[1],"%d",&startLevel);

    int maxLevel;
    sscanf(argv[2],"%d",&maxLevel);

    if( maxLevel < startLevel ){
      std::cout << "maxLevel >= startLevel not fulfilled." << std::endl;
    }

    // sequential version
    if (1 && helper.size()==1)
    {
#if HAVE_UG
      // make grid
      const int dim = 2;
      typedef Dune::UGGrid<dim> Grid;
      Dune::FieldVector<Grid::ctype, dim> ll(0.0);
      Dune::FieldVector<Grid::ctype, dim> ur(1.0);
      std::array<unsigned int, dim> elements;
      std::fill(elements.begin(), elements.end(), 1);

      std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
      grid->globalRefine(startLevel);

      typedef Grid::LeafGridView GV;
      GV gv = grid->leafGridView();
      adaptivity(*grid,gv,startLevel,maxLevel);
#else
      std::cout << "This example requires UG!" << std::endl;
#endif
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
