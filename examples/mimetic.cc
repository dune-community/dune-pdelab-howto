// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/float_cmp.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/mimeticfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridfunctionspace/intersectionindexset.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/poisson.hh>

#include "gridexamples.hh"

//===============================================================
// Problem setup and solution 
//===============================================================

// check DOFs in intersections
template<typename GV> 
void mimetictest (const GV& gv, std::string filename)
{
  // set up index set for intersections
  typedef Dune::PDELab::IntersectionIndexSet<GV> IIS;
  IIS iis(gv);

  // make finite element map
  typedef Dune::PDELab::MimeticLocalFiniteElementMap<IIS,double,double,GV::dimension> FEM;
  FEM fem(iis,Dune::GeometryType::cube);

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::NoConstraints,Dune::PDELab::ISTLVectorBackend<1>,
    Dune::PDELab::GridFunctionStaticSize<IIS> > GFS;
  GFS gfs(gv,fem,iis);

  // check result
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator; 
  typedef typename GV::IntersectionIterator IntersectionIterator;
  std::cout << std::endl;
  std::cout << "size of intersection index set = " << iis.size() << std::endl;
  for (ElementIterator it = gv.template begin<0>(); 
       it!=gv.template end<0>(); ++it)
    {
      //      std::cout << "element " << gv.indexSet().index(*it) << " has " << iis.size(*it) << " intersections" << std::endl;
      //   std::cout << "finite element map returns " << fem.find(*it).localCoefficients().size() << std::endl;

      IntersectionIterator endit = gv.iend(*it);
      int count=0;
      for (IntersectionIterator iit = gv.ibegin(*it); iit!=endit; ++iit)
        {
          //          std::cout << "  intersection " << count << " mapped to " << iis.subIndex(*it,count) << std::endl;
          if (iis.subIndex(*it,count)!=iis.index(*iit))
            std::cout << "error! go back to work on it!" << std::endl;
          count++;
        }
    }

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
  vtkwriter.write(filename,Dune::VTKOptions::ascii);
}


//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(1);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView(); 

      // solve problem
      mimetictest(gv,"mimetic_yasp_2d");
    }

#if HAVE_ALUGRID
    if (false)
    {
      // make grid 
//       ALUUnitSquare grid;
//       grid.globalRefine(0);

//       typedef ALUUnitSquare::Codim<0>::Partition<Dune::All_Partition>::LeafIterator Iterator;
//       typedef ALUUnitSquare::LeafIntersectionIterator IntersectionIterator;
//       typedef ALUUnitSquare::LeafGridView GV;
//       typedef ALUUnitSquare::ctype ctype;

//       // Do some random refinement. The result is a grid that may
//       // contain multiple hanging nodes per edge.
//       for(int i=0; i<4;++i){
//         Iterator it = grid.leafbegin<0,Dune::All_Partition>();
//         Iterator eit = grid.leafend<0,Dune::All_Partition>();

//         //        grid.mark(1,*(it));

//         for(;it!=eit;++it){
//           if((double)rand()/(double)RAND_MAX > 0.6)
//             grid.mark(1,*(it));
//         }
//         grid.preAdapt();
//         grid.adapt();
//         grid.postAdapt();
//       }

//       // get view
//       typedef ALUUnitSquare::LeafGridView GV;
//       const GV& gv=grid.leafView(); 
      
//       // solve problem
//       mimetictest(gv,"mimetic_ALU_2d");
    }

    // unit cube with hanging node refinement
    {
      // make grid 
      ALUCubeUnitSquare grid;
      grid.globalRefine(1);
      
      typedef ALUCubeUnitSquare::Codim<0>::Partition<Dune::All_Partition>::LeafIterator 
        Iterator;
      typedef ALUCubeUnitSquare::LeafIntersectionIterator IntersectionIterator;
      typedef ALUCubeUnitSquare::LeafGridView GV;
      typedef ALUCubeUnitSquare::ctype ctype;

      // get view
      const GV& gv=grid.leafView(); 
      const int dim = GV::dimension;
      
      // Do some random refinement. The result is a grid that may
      // contain multiple hanging nodes per edge.
      for(int i=0; i<4;++i){
        Iterator it = grid.leafbegin<0,Dune::All_Partition>();
        Iterator eit = grid.leafend<0,Dune::All_Partition>();

        for(;it!=eit;++it){
          if((double)rand()/(double)RAND_MAX > 0.6)
            grid.mark(1,*(it));
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();
      }

      // solve problem
      mimetictest(gv,"mimetic_ALUCUBE_3d");
    }

#endif

	// test passed
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
