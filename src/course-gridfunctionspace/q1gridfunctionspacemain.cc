/** \file \brief Set up a gridfunctionspace from local description */
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<dune/common/mpihelper.hh>
#include<dune/grid/yaspgrid.hh>
#include"q1gridfunctionspace.hh"
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // make grid
    Dune::FieldVector<double,2> L(1.0);
    Dune::FieldVector<int,2> N(2);
    Dune::FieldVector<bool,2> B(false);
    Dune::YaspGrid<2> grid(L,N,B,0);
    q1GridFunctionSpace(grid.leafView());

    return 0;
  }
  catch (Dune::Exception &e){ std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){ std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
