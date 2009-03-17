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
#endif
#include "basicunitcube.hh"

class YaspUnitSquare : public Dune::YaspGrid<2,2>
{
public:
  YaspUnitSquare () : Dune::YaspGrid<2,2>(Dune::FieldVector<double,2>(1.0),
					  Dune::FieldVector<int,2>(1),
					  Dune::FieldVector<bool,2>(false),0)
  {}
};

#if HAVE_ALUGRID
class ALUUnitSquare : public Dune::ALUSimplexGrid<2,2> 
{
public:
  ALUUnitSquare () : Dune::ALUSimplexGrid<2,2>("grids/2dsimplex.alu") {}
};

// class ALUReentrantCorner : public Dune::GridPtr<Dune::ALUSimplexGrid<2,2> >
// {
// public:
//   ALUReentrantCorner()
//     : Dune::GridPtr<Dune::ALUSimplexGrid<2,2> >("grids/2dreentrantcorner.dgf")
//   { }
// };
#endif //HAVE_ALUGRID


#if HAVE_ALBERTA
class AlbertaLDomain : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaLDomain () : Dune::AlbertaGrid<2,2>("grids/ldomain.al") {}
};

class AlbertaUnitSquare : public Dune::AlbertaGrid<2,2>
{
public:
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
class UGUnitSquare : public Dune::UGGrid<2>
{
public:
  UGUnitSquare (unsigned int heapSize=100) : Dune::UGGrid<2>(heapSize)
  {
	this->createBegin();
	Dune::FieldVector<double,2> pos;
	pos[0] = 0;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 1;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 0;  pos[1] = 1;
	this->insertVertex(pos);
	pos[0] = 1;  pos[1] = 1;
	this->insertVertex(pos);
	std::vector<unsigned int> cornerIDs(3);
	cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	this->createEnd();
  }
};

class UGUnitSquareQ : public Dune::UGGrid<2>
{
public:
  UGUnitSquareQ (unsigned int heapSize=100) : Dune::UGGrid<2>(heapSize)
  {
	this->createBegin();
	Dune::FieldVector<double,2> pos;
	pos[0] = 0;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 1;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 0;  pos[1] = 1;
	this->insertVertex(pos);
	pos[0] = 1;  pos[1] = 1;
	this->insertVertex(pos);
	std::vector<unsigned int> cornerIDs(4);
	cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;  cornerIDs[3] = 3;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);
	this->createEnd();
  }
};


class UGUnitTriangle : public Dune::UGGrid<2>
{
public:
  UGUnitTriangle (unsigned int heapSize=100) : Dune::UGGrid<2>(heapSize)
  {
	this->createBegin();
	Dune::FieldVector<double,2> pos;
	pos[0] = 0;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 1;  pos[1] = 0;
	this->insertVertex(pos);
	pos[0] = 0;  pos[1] = 1;
	this->insertVertex(pos);
	std::vector<unsigned int> cornerIDs(3);
	cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	this->createEnd();
  }
};

class UGLDomain : public Dune::UGGrid<2>
{
public:
  UGLDomain (unsigned int heapSize=100) : Dune::UGGrid<2>(heapSize)
  {
	this->createBegin();
	Dune::FieldVector<double,2> pos;
	pos[0] =-1.0;  pos[1] =-1.0; this->insertVertex(pos);
	pos[0] = 0.0;  pos[1] =-1.0; this->insertVertex(pos);
	pos[0] =-1.0;  pos[1] = 0.0; this->insertVertex(pos);
	pos[0] = 0.0;  pos[1] = 0.0; this->insertVertex(pos);
	pos[0] = 1.0;  pos[1] = 0.0; this->insertVertex(pos);
	pos[0] =-1.0;  pos[1] = 1.0; this->insertVertex(pos);
	pos[0] = 0.0;  pos[1] = 1.0; this->insertVertex(pos);
	pos[0] = 1.0;  pos[1] = 1.0; this->insertVertex(pos);
	std::vector<unsigned int> cornerIDs(3);
	cornerIDs[0] = 0;  cornerIDs[1] = 1;  cornerIDs[2] = 2;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 2;  cornerIDs[1] = 1;  cornerIDs[2] = 3;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 2;  cornerIDs[1] = 3;  cornerIDs[2] = 5;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 5;  cornerIDs[1] = 3;  cornerIDs[2] = 6;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 3;  cornerIDs[1] = 4;  cornerIDs[2] = 6;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	cornerIDs[0] = 6;  cornerIDs[1] = 4;  cornerIDs[2] = 7;
	this->insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);
	this->createEnd();
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
