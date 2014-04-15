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
#include <dune/common/parallel/mpihelper.hh>
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
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/stokesdg.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
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
    static const unsigned int vBlockSize = Dune::MonomImp::Size<dim,vOrder>::val;
    static const unsigned int pBlockSize = Dune::MonomImp::Size<dim,pOrder>::val;

    typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::no_blocking,vBlockSize> VVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::no_blocking,pBlockSize> PVectorBackend;
    typedef Dune::PDELab::ISTLVectorBackend<> VelocityVectorBackend;

#if 1
    // this creates a flat backend (i.e. blocksize == 1)
    typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::no_blocking> VectorBackend;
#else
    // this creates a backend with static blocks matching the size of the LFS
    typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::static_blocking> VectorBackend;
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

#if 0
    // enable for output of DOFIndices and mapped container indices
    {
      Dune::PDELab::LocalFunctionSpace<GFS> lfs(gfs);
      for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
        {compile
          lfs.bind(*it);
          for (int i = 0; i < lfs.size(); ++i)
            std::cout << lfs.dofIndex(i) << "  " << gfs.ordering().map_index(lfs.dofIndex(i)) << std::endl;
          std::cout << std::endl;
        }
    }
#endif

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
    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(75); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
    typedef Dune::PDELab::GridOperator
        <GFS,GFS,LocalDGOperator,MBE,RF,RF,RF,C,C> GOS;
    GOS gos(gfs,gfs,lop,mbe);

    std::cout << "=== grid operator space setup " <<  watch.elapsed() << " s" << std::endl;

    typedef typename GOS::Jacobian M;
    watch.reset();
    M m(gos);
    std::cout << m.patternStatistics() << std::endl;
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

    // // Create VTK Output
    // using namespace Dune::TypeTree;
    // typedef typename Dune::PDELab::GridFunctionSubSpace
    //   <GFS,TypePath<0> > velocitySubGFS;
    // velocitySubGFS velocitySubGfs(gfs);
    // typedef typename Dune::PDELab::GridFunctionSubSpace
    //   <GFS,TypePath<1> > pSubGFS;
    // pSubGFS pSubGfs(gfs);

    // typedef Dune::PDELab::VectorDiscreteGridFunction<velocitySubGFS,V> VDGF;
    // VDGF vdgf(velocitySubGfs,x);
    // typedef Dune::PDELab::DiscreteGridFunction<pSubGFS,V> PDGF;
    // PDGF pdgf(pSubGfs,x);

    //    #ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x);
    vtkwriter.write(filename,Dune::VTK::appendedraw);
    //    #endif
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
          const int dim = 2;
          Dune::FieldVector<double,dim> L(1.0);
          Dune::array<int,dim> N;
          N[0] = x; N[1] = y;
          std::bitset<dim> B(false);
          Dune::YaspGrid<dim> grid(L,N,B,0);

          // get view
            typedef Dune::YaspGrid<dim>::LeafGridView GV;
            const GV& gv=grid.leafGridView();

            // solve problem
            Dune::dinfo.push(false);
            stokes<GV,double,2,1>(gv,"dgstokes-2D-2-1", method);
        }


        // YaspGrid P2/P3 2D test
#if 1
        {
          // make grid
          const int dim = 3;
          Dune::FieldVector<double,dim> L(1.0);
          Dune::array<int,dim> N;
          N[0] = x; N[1] = y; N[2] = z;
          std::bitset<dim> B(false);
          Dune::YaspGrid<dim> grid(L,N,B,0);

          // get view
          typedef Dune::YaspGrid<dim>::LeafGridView GV;
          const GV& gv=grid.leafGridView();

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
