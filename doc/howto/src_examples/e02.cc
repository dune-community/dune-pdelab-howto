/** \file

    \brief Solve elliptic problem in constrained spaces with
    conforming finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include application heaeders
#include"e02_operator.hh"
#include"e02_parameter.hh"
#include"e02_P1.hh"

template <typename V>
void output(const V& v)
{
  for (int j=0; j<v.size(); ++j)
    std::cout << v[j] << "   ";
  std::cout << std::endl;
}

//===============================================================
// Problems
//===============================================================

// the Ljump is a building block of the eiffel tower!
void LEiffel(int level)
{
  // instanciate ug grid object
  typedef Dune::UGGrid<3> GridType;
  GridType grid(400);

  // vectors for boundary and material conditions
  std::vector<int> boundaryIndexToPhysicalEntity;
  std::vector<int> elementIndexToPhysicalEntity;

  // read a gmsh file
  std::string gridName = "./grids/L.msh";
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, gridName, boundaryIndexToPhysicalEntity,
      elementIndexToPhysicalEntity, true, false);

  // just for control
  output(elementIndexToPhysicalEntity);
  output(boundaryIndexToPhysicalEntity);

  // edit gridName
  gridName.erase(0, gridName.rfind("/")+1);
  gridName.erase(gridName.find(".", 0), gridName.length());

  // refine grid
  grid.globalRefine(level);

  // get a grid view
  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();

  // material and boundary conditions
  typedef LeiffelDiffusion<GV,double,std::vector<int> > M;
  M m(gv, elementIndexToPhysicalEntity);
  typedef LeiffelBCType<GV,std::vector<int> > B;
  B b(gv, boundaryIndexToPhysicalEntity);
  typedef LeiffelBCExtension<GV,double,std::vector<int> > G;
  G g(gv, boundaryIndexToPhysicalEntity);

  // Flux at boundaries
  typedef LeiffelFlux<GV,double,std::vector<int> > J;
  J j(gv, boundaryIndexToPhysicalEntity);

  // run simulation
  e02_P1(gv, m, b, g, j, gridName);
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

    // scan arguments
    if (argc!=2)
    {
      if(helper.rank()==0)
        std::cout << "usage: ./example02 <level>" << std::endl;
      return 1;
    }

    // refinement level
    int level;
    sscanf(argv[1],"%d",&level);

    // run simulations
    if (1 && helper.size()==1)
    {
      LEiffel(level);
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
