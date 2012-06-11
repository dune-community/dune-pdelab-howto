// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Stokes with DG method (stationary case).
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/entityblockedlocalordering.hh>
#include <dune/pdelab/gridfunctionspace/leaflocalordering.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/stokesdg.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include "sproblemA.hh"

#define USE_SUPER_LU
#define MAKE_VTK_OUTPUT

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename RF, int vOrder, int pOrder>
void stokes (const GV& gv, std::string filename, const std::string method)
{
    // <<<1>>> constants and types
    typedef typename GV::Grid::ctype DF;
    static const unsigned int dim = GV::dimension;
    Dune::Timer watch;
    std::cout << "=== Initialize" << std::endl;

    // <<<2>>> Make grid function space
    watch.reset();
    typedef Dune::PDELab::MonomLocalFiniteElementMap<DF,RF,dim,vOrder> vFEM;
    typedef Dune::PDELab::MonomLocalFiniteElementMap<DF,RF,dim,pOrder> pFEM;

    vFEM vFem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    pFEM pFem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    // DOFs per cell
#if 0
    static const unsigned int vBlockSize = Dune::MonomImp::Size<dim,vOrder>::val;
    static const unsigned int pBlockSize = Dune::MonomImp::Size<dim,pOrder>::val;
    static const unsigned int blockSize = dim * vBlockSize + pBlockSize;
    typedef Dune::PDELab::ISTLVectorBackend<
        Dune::PDELab::ISTLParameters::static_blocking, vBlockSize> VVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<
        Dune::PDELab::ISTLParameters::static_blocking, pBlockSize> PVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<
        Dune::PDELab::ISTLParameters::static_blocking, dim*vBlockSize> VelocityVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<
        Dune::PDELab::ISTLParameters::static_blocking, blockSize> VectorBackend;
#else
    typedef Dune::PDELab::ISTLVectorBackend<> VVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> PVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> VelocityVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> VectorBackend;
#endif
    // velocity
    Dune::dinfo << "--- v^dim" << std::endl;
    typedef Dune::PDELab::EntityBlockedOrderingTag VelocityOrderingTag;
    typedef Dune::PDELab::VectorGridFunctionSpace<
        GV,vFEM,dim,
        VelocityVectorBackend,
        VVectorBackend,
        Dune::PDELab::NoConstraints,
        VelocityOrderingTag
      > velocityGFS;
    velocityGFS velocityGfs(gv,vFem);
    velocityGfs.name("v");
    // p
    Dune::dinfo << "--- p" << std::endl;
    typedef Dune::PDELab::GridFunctionSpace<GV,pFEM,
        Dune::PDELab::NoConstraints, PVectorBackend> pGFS;
    pGFS pGfs(gv,pFem);
    pGfs.name("p");
    // GFS
    Dune::dinfo << "--- v^dim,p" << std::endl;
    typedef Dune::PDELab::EntityBlockedOrderingTag StokesOrderingTag;
    typedef Dune::PDELab::CompositeGridFunctionSpace<
        VectorBackend, StokesOrderingTag,
        velocityGFS, pGFS> GFS;
    GFS gfs(velocityGfs, pGfs);
    std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

    // <<<3>>> Make coefficient Vector and initialize it from a function
    typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
    V x(gfs);

    typedef B_A<GV> BType;
    BType b(gv);
    typedef F_A<GV,RF> FType;
    FType f(gv);
    typedef V_A<GV,RF> VType;
    VType v(gv);
    typedef P_A<GV,RF> PType;
    PType p(gv);

    // <<<4>>> Make grid Function operator
    watch.reset();
    const double mu = 1.0;
    typedef typename Dune::PDELab::DefaultInteriorPenalty<RF> PenaltyTerm;
    PenaltyTerm ip_term(method,mu);

    typedef Dune::PDELab::StokesDGParameters<GV,RF,FType,BType,VType,PType,PenaltyTerm>
      LocalDGOperatorParameters;
    LocalDGOperatorParameters lop_params(method, mu, f,b,v,p, ip_term);
    typedef Dune::PDELab::StokesDG<LocalDGOperatorParameters> LocalDGOperator;
    const int superintegration_order = 0;
    LocalDGOperator lop(lop_params,superintegration_order);

    typedef Dune::PDELab::EmptyTransformation C;
    typedef Dune::PDELab::GridOperator
        <GFS,GFS,LocalDGOperator, Dune::PDELab::ISTLMatrixBackend,RF,RF,RF,C,C> GOS;
    GOS gos(gfs,gfs,lop);

    std::cout << "=== grid operator space setup " <<  watch.elapsed() << " s" << std::endl;

    typedef typename GOS::Jacobian M;
    watch.reset();
    M m(gos);
    std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
    m = 0.0;
    watch.reset();
    gos.jacobian(x,m);
    std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;

    std::ofstream matrix("Matrix");
    Dune::printmatrix(matrix, m.base(), "M", "r", 6, 3);

    // evaluate residual w.r.t initial guess
    V r(gfs);
    r = 0.0;
    watch.reset();
    x = 0.0;
    gos.residual(x,r);
    std::cout << "=== residual evaluation " <<  watch.elapsed() << " s" << std::endl;

    bool verbose = true;

    typedef typename M::BaseT ISTLM;
    typedef typename V::BaseT ISTLV;
    #ifdef USE_SUPER_LU // use lu decomposition as solver
    #if HAVE_SUPERLU
    // make ISTL solver
    Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
    Dune::SuperLU<ISTLM> solver(m.base(), verbose?1:0);
    Dune::InverseOperatorResult stat;
    #else
    #error No superLU support, please install and configure it.
    #endif
    #else // Use iterative solver
    // make ISTL solver
    Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> opa(m.base());
    Dune::SeqILU0<ISTLM,ISTLV,ISTLV> ilu0(m.base(),1.0);
    typedef typename M::BaseT ISTLM;
    Dune::BiCGSTABSolver<ISTLV> solver(opa,ilu0,1E-10,20000, verbose?2:1);
    Dune::InverseOperatorResult stat;
    #endif

    // solve the jacobian system
    r *= -1.0; // need -residual
    x = r;
    solver.apply(x.base(),r.base(),stat);

    #ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    Dune::PDELab::add_solution_to_vtk_writer(
        vtkwriter, gfs, x);
    vtkwriter.write(filename,Dune::VTKOptions::binaryappended);
    #endif
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
    try{
        //Maybe initialize Mpi
        Dune::MPIHelper::instance(argc, argv);

        int x=2;
        int y=2;
        int z=2;
        std::string method = "nipg";
        if (argc > 1)
            x = atoi(argv[1]);
        if (argc > 2)
            y = atoi(argv[2]);
        if (argc > 3)
            z = atoi(argv[3]);
        if (argc > 4)
            method = argv[4];

        // YaspGrid P1/P2 2D test
        if(1){
            // make grid
            Dune::FieldVector<double,2> L(1.0);
            Dune::FieldVector<int,2> N(x); N[1] = y;
            Dune::FieldVector<bool,2> B(false);
            Dune::YaspGrid<2> grid(L,N,B,0);

            // get view
            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafView();

            // make finite element map
            typedef GV::Grid::ctype DF;
            // solve problem
            Dune::dinfo.push(false);
            stokes<GV,double,2,1>(gv,"dgstokes-2D-2-1", method);
        }


        // YaspGrid P2/P3 2D test
#if 0
        {
            // make grid
            Dune::FieldVector<double,3> L(1.0);
            Dune::FieldVector<int,3> N(x); N[1] = y; N[2] = z;
            Dune::FieldVector<bool,3> B(false);
            Dune::YaspGrid<3> grid(L,N,B,0);

            // get view
            typedef Dune::YaspGrid<3>::LeafGridView GV;
            const GV& gv=grid.leafView();

            // make finite element map
            typedef GV::Grid::ctype DF;
            // solve problem
            stokes<GV,double,2,1>(gv,"dgstokes-2D-3-2", method);
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
