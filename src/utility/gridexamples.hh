// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_GRIDEXAMPLES_HH
#define DUNE_PDELAB_GRIDEXAMPLES_HH

#include<dune/grid/yaspgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include<dune/grid/alugrid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#include "basicunitcube.hh"

// unit cubes/squares from test-parallel-ug
template <class Grid>
void insertVertices(int n, Dune::GridFactory<Grid> &factory, int coordIdx=0)
{
  Dune::FieldVector<double,Grid::dimension> pos;
  if (coordIdx == Grid::dimension) {
    factory.insertVertex(pos);
    return;
  }

  for (int i=0; i < n; ++i) {
    pos[coordIdx] = double(i)/(n-1);
    insertVertices(n, factory, coordIdx + 1);
  }
}

template <class GridType>
void insertElements(int n, Dune::GridFactory<GridType> &factory, unsigned coordIdx=0)
{
  const int dim = GridType::dimension;

  if (dim == 3) {
    for (int i=0; i<n-1; i++) {
      for (int j=0; j<n-1; j++) {
        for (int k=0; k<n-1; k++) {
          std::vector<unsigned int> v(8);
          v[0] = k*n*n + i*n + j;
          v[1] = k*n*n + i*n + j+1;
          v[2] = k*n*n + (i+1)*n + j;
          v[3] = k*n*n + (i+1)*n + j+1;

          v[4] = (k+1)*n*n + i*n + j;
          v[5] = (k+1)*n*n + i*n + j+1;
          v[6] = (k+1)*n*n + (i+1)*n + j;
          v[7] = (k+1)*n*n + (i+1)*n + j+1;

          factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), v);
        }
      }
    }
  }
  else if (dim == 2) {
    for (int i=0; i<n-1; i++) {
      for (int j=0; j<n-1; j++) {
        std::vector<unsigned int> v(4);
        v[0] = i*n + j;
        v[1] = i*n + j+1;
        v[2] = (i+1)*n + j;
        v[3] = (i+1)*n + j+1;

        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), v);
      }
    }
  }
}

template <class GridType>
void createGrid(GridType &grid, int n)
{
  typedef Dune::GridFactory<GridType> GridFactory;

  GridFactory factory(&grid);

  insertVertices<GridType>(n, factory);
  insertElements<GridType>(n, factory);

  factory.createGrid();
}

class YaspUnitSquare : public Dune::YaspGrid<2>
{
public:
  YaspUnitSquare () : Dune::YaspGrid<2>(Dune::FieldVector<double,2>(1.0),
                                        Dune::array<int,2>(Dune::fill_array<int,2>(1)),
                                        std::bitset<2>(false),
                                        0)
  {}
};

#if HAVE_ALUGRID
class ALUUnitSquare :
  public Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
{
public:
  ALUUnitSquare () :
    Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>("grids/2dsimplex.alu") {}
};

class ALUCubeUnitSquare : public Dune::ALUGrid<3,3,Dune::cube,Dune::conforming>
{
public:
  ALUCubeUnitSquare () :
    Dune::ALUGrid<3,3,Dune::cube,Dune::conforming>
    ("grids/3drefinedcube.alu") {}
};

// class ALUReentrantCorner : public Dune::GridPtr<Dune::ALUSimplexGrid<2,2> >
// {
// public:
//   ALUReentrantCorner()
//     : Dune::GridPtr<Dune::ALUSimplexGrid<2,2> >("grids/2dreentrantcorner.dgf")
//   { }
// };

template<int dim>
class ALUUnitCube : public BasicUnitCube< dim >
{
public:
  typedef Dune::ALUGrid<dim,dim,Dune::simplex,Dune::nonconforming> GridType;

private:
  GridType* grid_;

public:
  ALUUnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< dim >::insertVertices( factory );
    BasicUnitCube< dim >::insertSimplices( factory );
    grid_ = factory.createGrid();
  }

  ~ALUUnitCube ()
  {
    delete grid_;
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

#endif //HAVE_ALUGRID


#if HAVE_ALBERTA
class AlbertaLDomain : public Dune::AlbertaGrid<2,2>
{
public:
  typedef Dune::AlbertaGrid<2,2> Grid;
  AlbertaLDomain () : Dune::AlbertaGrid<2,2>("grids/ldomain.al") {}
};

class AlbertaUnitSquare : public Dune::AlbertaGrid<2,2>
{
public:
  typedef Dune::AlbertaGrid<2,2> Grid;
  AlbertaUnitSquare () : Dune::AlbertaGrid<2,2>("grids/2dgrid.al") {}
};

class AlbertaReentrantCorner : public Dune::GridPtr<Dune::AlbertaGrid<2,2> >
{
public:
  typedef Dune::AlbertaGrid<2,2> Grid;

  AlbertaReentrantCorner()
    : Dune::GridPtr<Dune::AlbertaGrid<2,2> >("grids/2dreentrantcorner.dgf")
  { }
};
#endif

#if HAVE_UG

struct UGSetHeapSize
{
  UGSetHeapSize (unsigned int heapSize)
  {
    Dune::UGGrid<2>::setDefaultHeapSize(heapSize);
  }
};

class UGUnitSquare : UGSetHeapSize, public Dune::UGGrid<2>
{
public:
  UGUnitSquare (unsigned int heapSize=100) : UGSetHeapSize(heapSize), Dune::UGGrid<2>()
  {
    Dune::GridFactory<Dune::UGGrid<2> > factory(this);
    Dune::FieldVector<double,2> pos;
    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;
    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;
    factory.insertVertex(pos);
    std::vector<unsigned int> cornerIDs(3);
    cornerIDs[0] = 1;  cornerIDs[1] = 3;  cornerIDs[2] = 0;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 2;  cornerIDs[1] = 0;  cornerIDs[2] = 3;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    factory.createGrid();
  }
};

class UGUnitSquareQ : UGSetHeapSize, public Dune::UGGrid<2>
{
public:
  UGUnitSquareQ (unsigned int heapSize=100) : UGSetHeapSize(heapSize), Dune::UGGrid<2>()
  {
    Dune::GridFactory<Dune::UGGrid<2> > factory(this);
    Dune::FieldVector<double,2> pos;
    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;
    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 1;
    factory.insertVertex(pos);
    std::vector<unsigned int> cornerIDs(4);
    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;  cornerIDs[3] = 3;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);
    factory.createGrid();
  }
};


class UGUnitTriangle : UGSetHeapSize, public Dune::UGGrid<2>
{
public:
  UGUnitTriangle (unsigned int heapSize=100) : UGSetHeapSize(heapSize), Dune::UGGrid<2>()
  {
    Dune::GridFactory<Dune::UGGrid<2> > factory(this);
    Dune::FieldVector<double,2> pos;
    pos[0] = 0;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 1;  pos[1] = 0;
    factory.insertVertex(pos);
    pos[0] = 0;  pos[1] = 1;
    factory.insertVertex(pos);
    std::vector<unsigned int> cornerIDs(3);
    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    factory.createGrid();
  }
};

class UGLDomain : UGSetHeapSize, public Dune::UGGrid<2>
{
public:
  UGLDomain (unsigned int heapSize=100) : UGSetHeapSize(heapSize), Dune::UGGrid<2>()
  {
    Dune::GridFactory<Dune::UGGrid<2> > factory(this);
    Dune::FieldVector<double,2> pos;
    pos[0] =-1.0;  pos[1] =-1.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] =-1.0; factory.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 1.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 1.0; factory.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 1.0; factory.insertVertex(pos);
    std::vector<unsigned int> cornerIDs(3);
    cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 2;  cornerIDs[1] = 3;  cornerIDs[2] = 5;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 5;  cornerIDs[1] = 3;  cornerIDs[2] = 6;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 3;  cornerIDs[1] = 4;  cornerIDs[2] = 6;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    cornerIDs[0] = 6;  cornerIDs[1] = 4;  cornerIDs[2] = 7;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
    factory.createGrid();
  }
};


class UGLDomainCubes : UGSetHeapSize, public Dune::UGGrid<2>
{
public:
  UGLDomainCubes (unsigned int heapSize=100) : UGSetHeapSize(heapSize), Dune::UGGrid<2>()
  {
    Dune::GridFactory<Dune::UGGrid<2> > factory(this);
    Dune::FieldVector<double,2> pos;
    pos[0] =-1.0;  pos[1] =-1.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] =-1.0; factory.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 0.0; factory.insertVertex(pos);
    pos[0] =-1.0;  pos[1] = 1.0; factory.insertVertex(pos);
    pos[0] = 0.0;  pos[1] = 1.0; factory.insertVertex(pos);
    pos[0] = 1.0;  pos[1] = 1.0; factory.insertVertex(pos);
    std::vector<unsigned int> cornerIDs(4);

    cornerIDs[0] = 0; cornerIDs[1] = 1; cornerIDs[2] = 2; cornerIDs[3] = 3;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);
    cornerIDs[0] = 2; cornerIDs[1] = 3; cornerIDs[2] = 5; cornerIDs[3] = 6;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);
    cornerIDs[0] = 3; cornerIDs[1] = 4; cornerIDs[2] = 6; cornerIDs[3] = 7;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);

    factory.createGrid();

    /* 5 - 6 - 7
       |   |   |
       2 - 3 - 4
       |   |
       0 - 1
    */
  }
};



template< int dim, int variant >
class UGUnitCube : public BasicUnitCube< dim >
{
public:
  typedef Dune::UGGrid<dim> GridType;

private:
  GridType* grid_;

public:
  UGUnitCube ()
  {
    Dune::GridFactory< GridType > factory;
    BasicUnitCube< dim >::insertVertices( factory );
    if( variant == 1 )
      BasicUnitCube< dim >::insertCubes( factory );
    else if( variant == 2 )
      BasicUnitCube< dim >::insertSimplices( factory );
    else
      DUNE_THROW( Dune::NotImplemented, "Variant "
                  << variant << " of UG unit cube not implemented." );
    grid_ = factory.createGrid();
  }

  ~UGUnitCube ()
  {
    delete grid_;
  }

  GridType &grid ()
  {
    return *grid_;
  }
};

#endif // HAVE_UG

#endif
