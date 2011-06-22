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
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/localoperator/stokesdg.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include "sproblemA.hh"

// #define USE_SUPER_LU
#define MAKE_VTK_OUTPUT

//===============================================================
// Problem setup and solution 
//===============================================================

// generate a P1 function and output it
template<typename GV, typename RF, int vOrder, int pOrder> 
void stokes (const GV& gv, std::string filename)
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
    static const unsigned int blockSize = dim * vBlockSize + pBlockSize;
    typedef Dune::PDELab::ISTLVectorBackend<blockSize> VectorBackend;
    // v
    Dune::dinfo << "--- v" << std::endl;
    typedef Dune::PDELab::GridFunctionSpace<GV,vFEM,
        Dune::PDELab::NoConstraints, VectorBackend,
        Dune::PDELab::SimpleGridFunctionStaticSize> vGFS;
    vGFS vGfs(gv,vFem);
    // velocity
    Dune::dinfo << "--- v^dim" << std::endl;
    typedef Dune::PDELab::PowerGridFunctionSpace<vGFS,dim,
        Dune::PDELab::ComponentBlockwiseOrderingTag<vBlockSize> > velocityGFS;
    velocityGFS velocityGfs(vGfs);
    // p
    Dune::dinfo << "--- p" << std::endl;
    typedef Dune::PDELab::GridFunctionSpace<GV,pFEM,
        Dune::PDELab::NoConstraints, VectorBackend,
        Dune::PDELab::SimpleGridFunctionStaticSize> pGFS;
    pGFS pGfs(gv,pFem);
    // GFS
    Dune::dinfo << "--- v^dim,p" << std::endl;
    typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::ComponentBlockwiseOrderingTag<dim * vBlockSize, pBlockSize>,
        velocityGFS, pGFS> GFS;
    GFS gfs(velocityGfs, pGfs);
    std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

#if 0
    // Iterater Grid
    typedef typename GFS::LocalFunctionSpace LFS;
    LFS lfs(gfs);
    typename GV::template Codim<0>::Iterator it = gv.template begin<0>();
    std::cout << "BlockSize " << blockSize << "\n";
    for (; it != gv.template end<0>(); ++it)
    {
        static int i = 0;
        lfs.bind(*it);
        std::cout << "Element " << i << " has " << lfs.size() << " DOFs\n Gobal Indices:";
        for (size_t j=0; j<lfs.size(); j++)
            std::cout << "    " << lfs.globalIndex(j);
        std::cout << "\n";

        typedef typename LFS::template Child<0>::Type LFSVelocity;
        const LFSVelocity& lfsvelocity = lfs.template getChild<0>();
        typedef typename LFSVelocity::template Child<0>::Type LFSv0;
        const LFSv0& lfsv0 = lfsvelocity.template getChild<0>();
        typedef typename LFSVelocity::template Child<1>::Type LFSv1;
        const LFSv1& lfsv1 = lfsvelocity.template getChild<1>();
        typedef typename LFS::template Child<1>::Type LFSP;
        const LFSP& lfsp = lfs.template getChild<1>();

        std::cout << "Element " << i << " has " << lfsvelocity.size() << " Velocity DOFs\n Gobal Indices:";
        for (size_t j=0; j<lfsvelocity.size(); j++)
            std::cout << "    " << lfsvelocity.globalIndex(j);
        std::cout << "\n Local Indices:";
        for (size_t j=0; j<lfsvelocity.size(); j++)
            std::cout << "    " << lfsvelocity.localIndex(j);
        std::cout << "\n";

        std::cout << "Element " << i << " has " << lfsv0.size() << " v0 DOFs\n Gobal Indices:";
        for (size_t j=0; j<lfsv0.size(); j++)
            std::cout << "    " << lfsv0.globalIndex(j);
        std::cout << "\n Local Indices:";
        for (size_t j=0; j<lfsv0.size(); j++)
            std::cout << "    " << lfsv0.localIndex(j);
        std::cout << "\n";
                
        std::cout << "Element " << i << " has " << lfsv1.size() << " v1 DOFs\n Gobal Indices:";
        for (size_t j=0; j<lfsv1.size(); j++)
            std::cout << "    " << lfsv1.globalIndex(j);
        std::cout << "\n Local Indices:";
        for (size_t j=0; j<lfsv1.size(); j++)
            std::cout << "    " << lfsv1.localIndex(j);
        std::cout << "\n";
                
        std::cout << "Element " << i << " has " << lfsp.size() << " pressure DOFs\n Gobal Indices:";
        for (size_t j=0; j<lfsp.size(); j++)
            std::cout << "    " << lfsp.globalIndex(j);
        std::cout << "\n Local Indices:";
        for (size_t j=0; j<lfsp.size(); j++)
            std::cout << "    " << lfsp.localIndex(j);
        std::cout << "\n";

        std::cout << "-------------------\n";
        
        i++;        
    }
    
    // Iterater Grid 2
    typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> velocitySubGFS;
    velocitySubGFS velocitySubGfs(gfs);
    typedef typename velocitySubGFS::LocalFunctionSpace velocitySubLFS;
    velocitySubLFS sublfsvelocity(velocitySubGfs);
    it = gv.template begin<0>();
    std::cout << "BlockSize " << blockSize << "\n";
    for (; it != gv.template end<0>(); ++it)
    {
        static int i = 0;
        sublfsvelocity.bind(*it);

        std::cout << "Element " << i << " has " << sublfsvelocity.size() << " Velocity DOFs\n Gobal Indices:";
        for (size_t j=0; j<sublfsvelocity.size(); j++)
            std::cout << "    " << sublfsvelocity.globalIndex(j);
        std::cout << "\n Local Indices:";
        for (size_t j=0; j<sublfsvelocity.size(); j++)
            std::cout << "    " << sublfsvelocity.localIndex(j);
        std::cout << "\n";

        std::cout << "-------------------\n";
        
        i++;        
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
    typedef Dune::PDELab::StokesDG<FType,BType,VType,PType> LocalDGOperator;
    const std::string method = "nipg";
    const double mu = 1.0;
    typename Dune::PDELab::DefaultInteriorPenalty<RF> ip_factor(method,mu);
    LocalDGOperator lop(method, ip_factor, mu, f,b,v,p);

    typedef Dune::PDELab::EmptyTransformation C;
    typedef Dune::PDELab::GridOperator
      <GFS,GFS,LocalDGOperator,typename VectorBackend::MatrixBackend,RF,RF,RF,C,C> GOS;
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

    for (size_t i = 0; i < r.base().N(); i++)
        for (size_t j = 0; j < r.base()[i].N(); j++)
            std::cout << i << "," << j << "\t" << r.base()[i][j] << "\n";
    bool verbose = true;
    
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
    Dune::BiCGSTABSolver<V> solver(opa,ilu0,1E-10,20000, verbose?2:1);
    Dune::InverseOperatorResult stat;
    #endif

    // solve the jacobian system
    r *= -1.0; // need -residual
    x = r;
    solver.apply(x,r,stat);
    //std::cout << x << std::endl;

    for (size_t i = 0; i < x.base().N(); i++)
        for (size_t j = 0; j < x.base()[i].N(); j++)
            std::cout << i << "," << j << "\t" << x.base()[i][j] << "\n";
    
    // make discrete function object
    // Important! We have to get the subspaces via GridFunctionSubSpace
    // the original (pre-composition) ones are not possible!
    typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> velocitySubGFS;
    velocitySubGFS velocitySubGfs(gfs);
    typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,1> pSubGFS;
    pSubGFS pSubGfs(gfs);
    typedef Dune::PDELab::VectorDiscreteGridFunction<velocitySubGFS,V> VDGF;
    VDGF vdgf(velocitySubGfs,x);
    typedef Dune::PDELab::DiscreteGridFunction<pSubGFS,V> PDGF;
    PDGF pdgf(pSubGfs,x);

    #ifdef MAKE_VTK_OUTPUT
    // output grid function with SubsamplingVTKWriter
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,2);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<PDGF>(pdgf,"p"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VDGF>(vdgf,"v"));
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VType>(v,"v0"));
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
        if (argc > 1)
            x = atoi(argv[1]);
        if (argc > 2)
            y = atoi(argv[2]);
        if (argc > 3)
            z = atoi(argv[3]);

        // YaspGrid P1/P2 2D test
        if(1){
            // make grid
            Dune::FieldVector<double,2> L(1.0);
            Dune::FieldVector<int,2> N(x); N[1] = y;
            Dune::FieldVector<bool,2> B(false);
            Dune::YaspGrid<2> grid(L,N,B,0);
//            grid.globalRefine(4);

            // get view
            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafView(); 

            // make finite element map
            typedef GV::Grid::ctype DF;
            // solve problem
            Dune::dinfo.push(false);
            stokes<GV,double,2,1>(gv,"dgstokes-2D-2-1");
        }

    
        // YaspGrid P2/P3 2D test
        if(0){
            // make grid
            Dune::FieldVector<double,3> L(1.0);
            Dune::FieldVector<int,3> N(x); N[1] = y; N[2] = z;
            Dune::FieldVector<bool,3> B(false);
            Dune::YaspGrid<3> grid(L,N,B,0);
            grid.globalRefine(0);

            // get view
            typedef Dune::YaspGrid<3>::LeafGridView GV;
            const GV& gv=grid.leafView(); 

            // make finite element map
            typedef GV::Grid::ctype DF;
            // solve problem
            stokes<GV,double,2,1>(gv,"dgstokes-2D-3-2");
        }

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
