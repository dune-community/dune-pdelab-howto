// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include<dune/common/mpihelper.hh>
#include<dune/grid/yaspgrid.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>

#include"q1localbasis.hh"
#include"q1localcoefficients.hh"
#include"q1localinterpolation.hh"
#include"q1localfiniteelement.hh"
#include"q1localfiniteelementmap.hh"
#include"q1analyticfunction.hh"
#include"q1constraints.hh"
#include"laplacedirichlet.hh"
#include"gridexamples.hh"

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    // Yasp 2D example
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2,2> grid(L,N,B,0);
      grid.globalRefine(5);

      // get view
      typedef Dune::YaspGrid<2,2>::LeafGridView GV;
      const GV& gv=grid.leafView(); 

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Q1LocalFiniteElementMap<DF,double> FEM;
      FEM fem;

      laplacedirichlet<GV,FEM,Q1Constraints>(gv,fem,3,"laplace_yasp_Q1_2d");
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2,2> grid(L,N,B,0);
      grid.globalRefine(2);

      // get view
      typedef Dune::YaspGrid<2,2>::LeafGridView GV;
      const GV& gv=grid.leafView(); 

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::Q22DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;
  
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,6,"laplace_yasp_Q2_2d");
    }

    // YaspGrid Q1 3D test
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(1);
      Dune::FieldVector<bool,3> B(false);
      Dune::YaspGrid<3,3> grid(L,N,B,0);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3,3>::LeafGridView GV;
      const GV& gv=grid.leafView(); 

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;
  
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,6,"laplace_yasp_Q1_3d");
    }

    // UG Pk 2D test
#if HAVE_UG
    {
      // make grid 
      UGUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef UGUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafView(); 
 
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
  
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,q,"laplace_UG_Pk_2d");
    }
#endif

    // UG P1 3D test
#if HAVE_UG
    {
      // make grid
      UGUnitCube<3,2> ugc; 
      ugc.grid().globalRefine(3);

      // get view
      typedef UGUnitCube<3,2>::GridType::LeafGridView GV;
      const GV& gv=ugc.grid().leafView(); 
 
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      typedef Dune::PDELab::P1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;
  
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,2,"laplace_UG_P1_3d");
    }
#endif

#if HAVE_ALBERTA
    {
      // make grid 
      AlbertaUnitSquare grid;
      grid.globalRefine(8);
      
      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafView(); 
      
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
      
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,q,"laplace_Alberta_Pk_2d");
    }
#endif

#if HAVE_ALUGRID
    {
      // make grid 
      ALUUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef ALUUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafView(); 
      
      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef double R;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);
      
      // solve problem
      laplacedirichlet<GV,FEM,Dune::PDELab::
        ConformingDirichletConstraints>(gv,fem,q,"laplace_ALU_Pk_2d");
    }
#endif

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
