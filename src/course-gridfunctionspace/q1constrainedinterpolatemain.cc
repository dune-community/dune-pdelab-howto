// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Initialize a grid function with constraints from a given function 
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include<dune/common/mpihelper.hh>
#include<dune/grid/yaspgrid.hh>

#include"q1localfiniteelementmap.hh"
#include"analyticfunction.hh"
#include"q1constraints.hh"
#include"q1constraintsparameters.hh"
#include"q1constrainedinterpolate.hh"

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // make grid
    Dune::FieldVector<double,2> L(1.0);
    Dune::FieldVector<int,2> N(1);
    Dune::FieldVector<bool,2> B(false);
    Dune::YaspGrid<2> grid(L,N,B,0);
    grid.globalRefine(5);
    q1interpolate(grid.leafView());
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
