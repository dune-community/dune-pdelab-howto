// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Taylor-Hood gridfunctionspace example 
*/
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/yaspgrid.hh>

#include"analyticfunction.hh"
#include"thvelocity.hh"
#include"thinterpolate.hh"

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // make grid
    Dune::FieldVector<double,2> L(1.0);
    Dune::array<int,2> N(Dune::fill_array<int,2>(1));
    std::bitset<2> periodic(false);
    Dune::YaspGrid<2> grid(L,N,periodic,0);
    grid.globalRefine(5);
    thinterpolate(grid.leafGridView());
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
