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
#include<dune/common/parallel/mpihelper.hh>
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
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/localoperator/diffusiondg.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include"../utility/gridexamples.hh"

// Select Problem
#include"problemA.hh"  // exp(-norm(x,y))
#include"problemB.hh"  // Like problem A but corners have small parts with Neumann boundary
#include"problemC.hh"  // Constant flow with checkerboard changing permeability, kind of ground water problem
#include"problemD.hh"  // Constant flow with randomly changig permeability
#include"problemE.hh"  // Constant flow with constant permeability with is 1E-6 in x direction
#include"problemF.hh"  // Constant flow with constant permeability with is 1E-6 in x direction

// Define most changed values

#define ProblemC

static const int monom_max_order = 5;
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
  typedef double Real;
    typedef double RF;
    Dune::Timer watch;



    // choose problem parameters

#ifdef ProblemA
    typedef K_A<GV,RF> KType;
    typedef F_A<GV,RF> FType;
    typedef BCTypeParam_A BType;
    typedef G_A<GV,RF> GType;
    typedef J_A<GV,RF> JType;
#endif

#ifdef ProblemB
    typedef K_B<GV,RF> KType;
    typedef F_B<GV,RF> FType;
    typedef BCTypeParam_B BType;
    typedef G_B<GV,RF> GType;
    typedef J_B<GV,RF> JType;
#endif

#ifdef ProblemC
    typedef K_C<GV,RF> KType;
    typedef F_C<GV,RF> FType;
    typedef BCTypeParam_C BType;
    typedef G_C<GV,RF> GType;
    typedef J_C<GV,RF> JType;
#endif

    KType k(gv);
    FType f(gv);
    BType bctype;
    GType g(gv);
    JType j(gv);

    // make grid function space
    typedef Dune::PDELab::NoConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking,BLOCK_SIZE> VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

    watch.reset();
    GFS gfs(gv,fem);
    if (verbose)
    {
        std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;
    }

    // make local operator
    typedef Dune::PDELab::DiffusionDG<KType,FType,BType,GType,JType> LOP;
    LOP lop(k,f,bctype,g,j,DG_METHOD);

    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(27); // Number of diagonals is depending on dim, order and geometry.
    typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    CC cc;

    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,Real,Real,Real,CC,CC> GO;
    GO go(gfs,cc,gfs,cc,lop,mbe);

    // make a vector of degree of freedom vectors
    typedef typename GO::Traits::Domain V;
    V solution(gfs,0.0);

    //Warning: Iterative solver does not work. Please use SuperLU here!
    //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
    //LS ls(10000,2);
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls(1);
    typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,V> SLP;
    SLP slp(go,ls,solution,1e-12);
    slp.apply();

    // make discrete function object
    typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
    DGF dgf(gfs,solution);

#ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    //KType<GV,RF> k2(gv);
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"u"));
    //vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<KType>(k,"k"));
    vtkwriter.addCellData(new Dune::PDELab::VTKFiniteElementMapAdapter<GV,FEM>(fem,"fem order"));
    vtkwriter.write(filename,Dune::VTK::appendedraw);
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

    try
    {
        // 2D
        {
            // make grid
            const int dim=2;
            Dune::FieldVector<double,dim> L(1.0);  // L[0]=2.0; L[1]=1.0;
            Dune::array<int,dim> N(Dune::fill_array<int,dim>(16));
            std::bitset<dim> P(false);
            typedef Dune::YaspGrid<dim> Grid;
            Grid grid(L,N,P,0);
#ifdef REFINE_STEPWISE
            for (int i = 0; i <= GRID_REFINE; ++i)
            {
                solve_dg(grid.leafGridView(), false);
                grid.globalRefine(1);
            }
#else
            // grid.globalRefine(GRID_REFINE);

            // instantiate finite element maps
            typedef Dune::SingleCodimSingleGeomTypeMapper<Grid::LeafGridView, 0> CellMapper;
            CellMapper cellmapper(grid.leafGridView());
            typedef Dune::PDELab::VariableMonomLocalFiniteElementMap<
                CellMapper,double,double,2,monom_max_order> FEM;
            FEM fem(cellmapper, 2); // works only for cubes

            // set polynomial order per element
            unsigned int range = std::ceil( double(cellmapper.size()) / (monom_max_order-1) );
            for (Grid::LeafGridView::Codim<0>::Iterator it = grid.leafGridView().begin<0>(), end = grid.leafGridView().end<0>();
                 it != end; ++it)
            {
                fem.setOrder(*it, 2 + cellmapper.map(*it) / range);
            }

            // solve problem :)
            solve_dg(grid.leafGridView(),fem,"DG_Yasp_2d",true);
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
