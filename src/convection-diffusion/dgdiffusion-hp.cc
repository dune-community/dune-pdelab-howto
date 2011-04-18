// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Discontinuous Galerkin with variable polynomial degree
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
#include<dune/grid/common/scsgmapper.hh>

#include<dune/pdelab/finiteelementmap/variablemonomfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/diffusiondg.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

#include"../utility/gridexamples.hh"

// Select Problem
#include"problemA.hh"  // exp(-norm(x,y))
#include"problemB.hh"  // Like problem A but corners have small parts with Neumann boundary
#include"problemC.hh"  // Constant flow with checkerboard changing permeability, kind of ground water problem
#include"problemD.hh"  // Constant flow with randomly changig permeability
#include"problemE.hh"  // Constant flow with constant permeability with is 1E-6 in x direction
#include"problemF.hh"  // Constant flow with constant permeability with is 1E-6 in x direction

// Define most changed values
static const int monom_max_order = 10;
#define BLOCK_SIZE 1
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
    typedef typename GV::Grid Grid;
    typedef typename GV::Grid::ctype DF;
    typedef double RF;
    Dune::Timer watch;

    // make function space
    typedef Dune::PDELab::ISTLVectorBackend<BLOCK_SIZE> VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
       Dune::PDELab::NoConstraints,
       VBE> GFS;
    watch.reset();
    GFS gfs(gv,fem);
    if (verbose)
    {
        std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;
    }

    // make coefficient Vector and initialize it from a function
    typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
    
    V x0(gfs);
    x0 = 0.0;
    V v(gfs,0.0);
    typedef K_C<GV,RF> KType;
    KType k(gv);
    typedef F_C<GV,RF> FType;
    FType f(gv);

    typedef BCTypeParam_C BCType;
    BCType bctype;

    typedef G_C<GV,RF> GType;
    GType g(gv);
    typedef J_C<GV,RF> JType;
    JType j(gv);
    Dune::PDELab::interpolate(g,gfs,x0);

    // make grid function operator
    typedef Dune::PDELab::DiffusionDG<KType,FType,BCType,GType,JType> LOP;
    LOP la(k,f,bctype,g,j,DG_METHOD);
    typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,
       Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation,
       VBE::MatrixBackend> GOS;
    GOS gos(gfs,gfs,la);

    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(1);
    typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,V> SLP;
    SLP slp(gos,v,ls,1e-12);
    slp.apply();

    // make discrete function object
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF dgf(gfs,v);

#ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    k_C<GV,RF> k2(gv);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,5);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"u"));
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter< k_C<GV,RF> >(k2,"k"));
    vtkwriter.addCellData(new Dune::PDELab::VTKFiniteElementMapAdapter<Grid,FEM>(fem,"fem order"));
    vtkwriter.write(filename,Dune::VTKOptions::binaryappended);
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
        {
            // make grid
            Dune::FieldVector<double,2> L(1.0);  // L[0]=2.0; L[1]=1.0;
            Dune::FieldVector<int,2> N(32);       // N[0]=2; N[1]=2;
            Dune::FieldVector<bool,2> B(false);
            typedef Dune::YaspGrid<2> Grid;
            Grid grid(L,N,B,0);
            #ifdef REFINE_STEPWISE
            for (int i = 0; i <= GRID_REFINE; ++i)
            {
                solve_dg(grid.leafView(), false);
                grid.globalRefine(1);
            }
            #else
            // grid.globalRefine(GRID_REFINE);

            // instantiate finite element maps
            typedef Dune::SingleCodimSingleGeomTypeMapper<Grid::LeafGridView, 0> CellMapper;
            CellMapper cellmapper(grid.leafView());
            typedef Dune::PDELab::VariableMonomLocalFiniteElementMap<
                CellMapper,double,double,2,monom_max_order> FEM;
            FEM fem(cellmapper, 2); // works only for cubes

            // set polynomial order per element
            unsigned int range = std::ceil( double(cellmapper.size()) / (monom_max_order-1) );
            for (Grid::LeafGridView::Codim<0>::Iterator it = grid.leafView().begin<0>(), end = grid.leafView().end<0>();
                 it != end; ++it)
            {
                fem.setOrder(*it, 2 + cellmapper.map(*it) / range);
            }
            
            // solve problem :)
            solve_dg(grid.leafView(),fem,"DG_Yasp_2d",true);
            #endif
        }

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
