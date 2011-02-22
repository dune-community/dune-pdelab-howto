// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Discontinuous Galerkin method 
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>
#include<dune/istl/superlu.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/diffusiondg.hh>

#include"../utility/gridexamples.hh"

// Select Problem
#include"problemA.hh"  // exp(-norm(x,y))
#include"problemB.hh"  // Like problem A but corners have small parts with Neumann boundary
#include"problemC.hh"  // Constant flow with checkerboard changing permeability, kind of ground water problem
#include"problemD.hh"  // Constant flow with randomly changig permeability
#include"problemE.hh"  // Constant flow with constant permeability with is 1E-6 in x direction
#include"problemF.hh"  // Constant flow with constant permeability with is 1E-6 in x direction

// Define most changed values
#define MONOM_BASIS_ORDER 2
#define BLOCK_SIZE 6
#define GRID_REFINE 2
#define DG_METHOD 0  // OBB: 0, NIPG: 1, SIPG: 2
#define MAKE_VTK_OUTPUT
//#define CALCULATE_L2_ERROR
//#define CALCULATE_ABSOLUTE_ERROR
#define USE_SUPER_LU
//#define REFINE_STEPWISE

// Solve the given problem with OBB, NIPG or SIPG
template<class GV, class FEM>
void solve_dg (const GV& gv, const FEM& fem, std::string filename, const bool verbose)
{
    typedef typename GV::Grid::ctype DF;
    typedef double RF;
    Dune::Timer watch;

    // make function space
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
       Dune::PDELab::NoConstraints,
       Dune::PDELab::ISTLVectorBackend<BLOCK_SIZE> > GFS;
    watch.reset();
    GFS gfs(gv,fem);
    if (verbose)
    {
        std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;
    }

    // make coefficient Vector and initialize it from a function
    typedef typename GFS::template VectorContainer<RF>::Type V;
    V x0(gfs);
    x0 = 0.0;
    typedef K_C<GV,RF> KType;
    KType k(gv);
    typedef F_C<GV,RF> FType;
    FType f(gv);
    typedef B_C<GV> BType;
    BType b(gv);
    typedef G_C<GV,RF> GType;
    GType g(gv);
    typedef J_C<GV,RF> JType;
    JType j(gv);
    Dune::PDELab::interpolate(g,gfs,x0);

    // make grid function operator
    typedef Dune::PDELab::DiffusionDG<KType,FType,BType,GType,JType> LOP;
    LOP la(k,f,b,g,j,DG_METHOD);
    typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,
       Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,
       Dune::PDELab::ISTLBCRSMatrixBackend<BLOCK_SIZE,BLOCK_SIZE> > GOS;
    GOS gos(gfs,gfs,la);

    // represent operator as a matrix
    typedef typename GOS::template MatrixContainer<RF>::Type M;
    watch.reset();
    M m(gos);
    if (verbose)
    {
        std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
    }
    m = 0.0;
    watch.reset();
    gos.jacobian(x0,m);
    if (verbose)
    {
        std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
    }
    //Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

    // evaluate residual w.r.t initial guess
    V r(gfs);
    r = 0.0;
    watch.reset();
    x0 = 0.0;
    gos.residual(x0,r);
    //std::cout << "Residuenvektor" << std::endl << r << std::endl;
    if (verbose)
    {
        std::cout << "=== residual evaluation " <<  watch.elapsed() << " s" << std::endl;
    }

    #ifdef USE_SUPER_LU // use lu decomposition as solver
    #if HAVE_SUPERLU
    // make ISTL solver
    Dune::MatrixAdapter<M,V,V> opa(m);
    typedef typename M::BaseT ISTLM;
    Dune::SuperLU<ISTLM> solver(m, verbose?1:0);
    Dune::InverseOperatorResult stat;
    #else
    #error No superLU support, please install and configure it.
    #endif
    #else // Use iterative solver
    // make ISTL solver
    Dune::MatrixAdapter<M,V,V> opa(m);
    Dune::SeqILU0<M,V,V> ilu0(m,1.0);
    typedef typename M::BaseT ISTLM;
    Dune::BiCGSTABSolver<V> solver(opa,ilu0,1E-10,20000,1);//verbose?1:0);
    Dune::InverseOperatorResult stat;
    #endif

    // solve the jacobian system
    r *= -1.0; // need -residual
    V x(gfs,0.0);
    solver.apply(x,r,stat);
    x += x0;
    //std::cout << x << std::endl;

    // make discrete function object
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF dgf(gfs,x);

    #ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"u"));
    vtkwriter.write(filename,Dune::VTKOptions::ascii);
    #endif
}

int main(int argc, char** argv)
{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
		if(helper.rank()==0)
		  std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    
    std::string problem="C";

    try
    {
        // 2D
      if (true)
        {
            // make grid
            Dune::FieldVector<double,2> L(1.0);  // L[0]=2.0; L[1]=1.0;
            Dune::FieldVector<int,2> N(64);       // N[0]=2; N[1]=2;
            Dune::FieldVector<bool,2> B(false);
            Dune::YaspGrid<2> grid(L,N,B,0);
            #ifdef REFINE_STEPWISE
            for (int i = 0; i <= GRID_REFINE; ++i)
            {
                solve_dg(grid.leafView(), false);
                grid.globalRefine(1);
            }
            #else
            grid.globalRefine(GRID_REFINE);

            // instantiate finite element maps
            typedef Dune::PDELab::MonomLocalFiniteElementMap<double,double,2,MONOM_BASIS_ORDER> FEM;
            FEM fem(Dune::GeometryType(Dune::GeometryType::cube,2)); // works only for cubes

            // solve problem :)
            solve_dg(grid.leafView(),fem,"DG_Yasp_2d",true);
            #endif
        }

#if HAVE_ALBERTA
      if (true)
        {
          typedef AlbertaUnitSquare GridType;
          GridType grid;
          grid.globalRefine(10);
          
          // get view
          typedef GridType::LeafGridView GV;
          const GV& gv=grid.leafView(); 
 
          // instantiate finite element maps
          typedef Dune::PDELab::MonomLocalFiniteElementMap<double,double,2,MONOM_BASIS_ORDER> FEM;
          FEM fem(Dune::GeometryType(Dune::GeometryType::simplex,2)); // works only for cubes

          solve_dg(gv,fem,"DG_Alberta_2d",true);
         
        }
#endif


    }
    catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
