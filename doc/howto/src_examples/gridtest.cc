// gridtest.cc: Read (sequential) gmsh files.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<string>

#include<dune/common/exceptions.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>      // New: Dune::GmshReader
#if HAVE_UG                                    // New: Use UG here
#include<dune/grid/uggrid.hh>
#endif

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  // check arguments
  if (argc!=3)
  {
    std::cout << "usage: " << argv[0] << "<gridName> <level>" << std::endl;
    return 1;
  }

  // refinement level
  int level = 0;
  sscanf(argv[2], "%d", &level);

  // instanciate ug grid object (xxx MB heap)
  typedef Dune::UGGrid<3> GridType;
  GridType grid(400);

  // read a gmsh file
  std::string gridName = argv[1];
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, gridName);

  // edit gridName
  gridName.erase(0, gridName.rfind("/")+1);
  gridName.erase(gridName.find(".",0), gridName.length());

  // refine grid
  grid.globalRefine(level);

  // get a grid view
  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();

  // plot celldata
  std::vector<int> a(gv.size(0), 1);

  //  output
  Dune::VTKWriter<GV> vtkwriter(gv);
  vtkwriter.addCellData(a, "celldata");
  vtkwriter.write(gridName.c_str(), Dune::VTKOptions::ascii);

  // done
  return 0;
}
