// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file 
    \brief Mimetic finite difference example (solves Problem A-F)
*/
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
#include<dune/pdelab/finiteelementmap/mfdconstraints.hh>
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
#include<dune/pdelab/localoperator/diffusionmfd.hh>
#include<dune/pdelab/localoperator/poisson.hh>

#include "../utility/gridexamples.hh"
#include "problemA.hh"
#include "problemB.hh"
#include "problemC.hh"
#include "problemD.hh"
#include "problemE.hh"
#include "problemF.hh"

template<typename K0, typename A0, typename F, typename B, typename J, typename G>
class DiffusionData
{
    typedef typename K0::Traits::GridViewType GV;
public:
    typedef GV GridView;
    static const unsigned int dimension = GV::dimension;
    typedef typename K0::Traits::DomainFieldType ctype;
    typedef typename K0::Traits::RangeFieldType rtype;

    typedef K0 KType;
    typedef A0 A0Type;
    typedef F FType;
    typedef B BType;
    typedef J JType;
    typedef G GType;

public:
    enum BCType { bcDirichlet, bcNeumann };

    DiffusionData(const GV& gv_)
        : kfunc(gv_), a0func(gv_), ffunc(gv_), bfunc(gv_), jfunc(gv_), gfunc(gv_), gv(gv_)
    {}

    template<class Entity>
    typename K0::Traits::RangeType K(const Entity& ent, const typename K0::Traits::DomainType& x) const
    {
        typename K0::Traits::RangeType y;
        kfunc.evaluate(ent, x, y);
        return y;
    }

    template<class Entity>
    typename A0::Traits::RangeType a_0(const Entity& ent, const typename A0::Traits::DomainType& x) const
    {
        typename A0::Traits::RangeType y;
        a0func.evaluate(ent, x, y);
        return y;
    }

    template<class Entity>
    typename F::Traits::RangeType f(const Entity& ent, const typename F::Traits::DomainType& x) const
    {
        typename F::Traits::RangeType y;
        ffunc.evaluate(ent, x, y);
        return y;
    }

    template<class Intersection>
    BCType bcType(const Intersection& is, const typename B::Traits::DomainType& x) const
    {
        typename B::Traits::RangeType y;
        bfunc.evaluate(is, x, y);
        if (y > 0)
            return bcDirichlet;
        else
            return bcNeumann;
    }

    template<class Intersection>
    typename J::Traits::RangeType j(const Intersection& is, const typename J::Traits::DomainType& x) const
    {
        typename J::Traits::RangeType y;
        jfunc.evaluate(is, x, y);
        return y;
    }

    const GV& gridview() const
    {
        return gv;
    }

    K0 kfunc;
    A0 a0func;
    F ffunc;
    B bfunc;
    J jfunc;
    G gfunc;

private:
    const GV& gv;
};

template<typename GV>
class Dummy : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
    BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,Dummy<GV> >
{
    const GV& gv;

public:
    typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
    typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,Dummy<GV> > BaseT;

    Dummy (const GV& gv_) : gv(gv_) {}

    template<typename I>
    inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                          const typename Traits::DomainType& x,
                          typename Traits::RangeType& y) const
    {}

    //! get a reference to the GridView
    inline const GV& getGridView ()
    {
        return gv;
    }
};

template<class B, class G, class GFS, class V>
void mimeticDirichletBoundaryConditions(const B& b, const G& g, const GFS& gfs, V& x)
{
    typedef typename GFS::Traits::GridViewType GV;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::Intersection Intersection;

    static const unsigned int dimIntersection = B::Traits::dimDomain;
    typedef typename B::Traits::DomainFieldType ctype;

    const GV& gridview = gfs.gridview();

    // make local function space
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,1> FaceSpace;
    FaceSpace fs(gfs);
    typedef typename FaceSpace::LocalFunctionSpace FaceUnknowns;
    FaceUnknowns face_space(fs);

    ElementIterator itend = gridview.template end<0>();
    for (ElementIterator it = gridview.template begin<0>(); it != itend; ++it)
    {
        // bind local function space to element
        face_space.bind(*it);

        unsigned int face = 0;
        IntersectionIterator isend = gridview.iend(*it);
        for (IntersectionIterator is = gridview.ibegin(*it); is != isend; ++is, ++face)
        {
            Dune::GeometryType gt = is->type();
            typename B::Traits::DomainType center
                = Dune::GenericReferenceElements<ctype,dimIntersection>::general(gt).position(0,0);
            typename B::Traits::RangeType bctype;
            b.evaluate(Dune::PDELab::IntersectionGeometry<Intersection>(*is, face), center, bctype);
            if (bctype > 0)
            {
                typename G::Traits::DomainType local_face_center
                    = is->geometryInInside().global(center);
                g.evaluate(*it, local_face_center, x[face_space.globalIndex(face)]);
            }
        }
    }
}

//===============================================================
// Problem setup and solution
//===============================================================

template<typename Data>
void mimetictest(Data& data, std::string filename)
{
    typedef typename Data::GridView GV;
    const GV& gv = data.gridview();

    const int dim = GV::dimension;

    // set up index set for intersections
    typedef Dune::PDELab::IntersectionIndexSet<GV> IIS;
    IIS iis(gv);

    // make finite element maps
    typedef Dune::PDELab::P0LocalFiniteElementMap<double,double,dim> CellFEM;
    CellFEM cell_fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    typedef Dune::PDELab::MimeticLocalFiniteElementMap<IIS,double,double,dim> FaceFEM;
    FaceFEM face_fem(iis, Dune::GeometryType::cube);

    // make function spaces
    typedef Dune::PDELab::GridFunctionSpace<GV,CellFEM,
        Dune::PDELab::NoConstraints,Dune::PDELab::ISTLVectorBackend<1> > CellGFS;
    CellGFS cell_gfs(gv, cell_fem);
    typedef Dune::PDELab::GridFunctionSpace<GV,FaceFEM,
        Dune::PDELab::MimeticConstraints,Dune::PDELab::ISTLVectorBackend<1>,
        Dune::PDELab::GridFunctionStaticSize<IIS> > FaceGFS;
    FaceGFS face_gfs(gv, face_fem, iis);
    typedef Dune::PDELab::CompositeGridFunctionSpace
        <Dune::PDELab::GridFunctionSpaceLexicographicMapper,CellGFS,FaceGFS> GFS;
    GFS gfs(cell_gfs, face_gfs);

    // construct a composite boundary condition type function
    typedef Dummy<GV> DType;
    DType d(gv);
    typedef typename Data::BType BType;
    typedef Dune::PDELab::CompositeGridFunction<DType,BType> BCT;
    BCT bct(d, data.bfunc);

    // constraints
    typedef typename GFS::template ConstraintsContainer<double>::Type T;
    T t;                               // container for transformation
    Dune::PDELab::constraints(bct,gfs,t); // fill container

    // make coefficent vector
    typedef typename GFS::template VectorContainer<double>::Type V;
    V x(gfs);

    // set Dirichlet boundary conditions
    mimeticDirichletBoundaryConditions(data.bfunc, data.gfunc, gfs, x);

    typedef Dune::PDELab::DiffusionMFD<Data> LOP;
    LOP lop(data);
    typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
        LOP,T,T,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
    GOS gos(gfs,t,gfs,t,lop);

    // evaluate residual w.r.t initial guess
    V r(gfs);
    r = 0.0;
    gos.residual(x,r);

    typedef typename GOS::template MatrixContainer<double>::Type M;
    M m(gos);
    m = 0.0;
    gos.jacobian(x,m);

    // make ISTL solver
    Dune::MatrixAdapter<M,V,V> opa(m);
    Dune::SeqSSOR<M,V,V> ssor(m,1,1.0);
    Dune::CGSolver<V> solver(opa,ssor,1E-10,5000,2);
    Dune::InverseOperatorResult stat;

    // solve the jacobian system
    V z(gfs,0.0);
    solver.apply(z,r,stat);
    x -= z;

    // make discrete function object
    typedef Dune::PDELab::GridFunctionSubSpace<GFS,0> PSUB;
    PSUB psub(gfs);                   // pressure subspace
    typedef Dune::PDELab::DiscreteGridFunction<PSUB,V> DGF;
    DGF dgf(psub, x);

    // output grid function with VTKWriter
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
    vtkwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
    vtkwriter.write(filename,Dune::VTKOptions::ascii);
}


//===============================================================
// Main program with grid setup
//===============================================================

template<class GV>
struct Data
{
    typedef DiffusionData<K_A<GV,double>, A0_A<GV,double>, F_A<GV,double>, B_A<GV>, J_A<GV,double>, G_A<GV,double> > A;
    typedef DiffusionData<K_B<GV,double>, A0_B<GV,double>, F_B<GV,double>, B_B<GV>, J_B<GV,double>, G_B<GV,double> > B;
    typedef DiffusionData<K_C<GV,double>, A0_C<GV,double>, F_C<GV,double>, B_C<GV>, J_C<GV,double>, G_C<GV,double> > C;
    typedef DiffusionData<K_D<GV,double>, A0_D<GV,double>, F_D<GV,double>, B_D<GV>, J_D<GV,double>, G_D<GV,double> > D;
    typedef DiffusionData<K_E<GV,double>, A0_E<GV,double>, F_E<GV,double>, B_E<GV>, J_E<GV,double>, G_E<GV,double> > E;
    typedef DiffusionData<K_F<GV,double>, A0_F<GV,double>, F_F<GV,double>, B_F<GV>, J_F<GV,double>, G_F<GV,double> > F;
};

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
            grid.globalRefine(5);

            // get view
            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafView();

            // choose problem data
            Data<GV>::C data(gv);

            // solve problem
            mimetictest(data,"mimetic_yasp_2d");
        }

#if HAVE_UG

        {
            // make grid
            typedef UGUnitSquareQ Grid;
            Grid grid;
            grid.setClosureType(Grid::NONE);
            grid.globalRefine(5);

            typedef Grid::Codim<0>::LeafIterator Iterator;
            typedef Iterator::Entity::Geometry Geometry;

            for(int i = 0; i < 1; ++i)
            {
                Iterator itend = grid.leafend<0>();
                for (Iterator it = grid.leafbegin<0>(); it != itend; ++it)
                {
                    const Geometry& geo = it->geometry();
                    for (int j = 0; j < geo.corners(); ++j)
                    {
                        Dune::FieldVector<Geometry::ctype,2> x = geo.corner(j);
                        double dummy;
                        if (std::modf(x[0]*8.0+1e-10, &dummy) < 2e-10 &&
                            std::modf(x[1]*8.0+1e-10, &dummy) < 2e-10)
                        {
                            grid.mark(1, *it);
                            break;
                        }
                    }
                }
                grid.preAdapt();
                grid.adapt();
                grid.postAdapt();
            }

            // get view
            typedef Grid::LeafGridView GV;
            const GV& gv=grid.leafView();

            Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::nonconforming);
            vtkwriter.write("ug",Dune::VTKOptions::ascii);

            // choose problem data
            Data<GV>::C data(gv);

            // solve problem
            mimetictest(data,"mimetic_ug_2d");
        }

#endif

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
//             ALUCubeUnitSquare grid;
//             grid.globalRefine(1);

//             typedef ALUCubeUnitSquare::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
//                 Iterator;
//             typedef ALUCubeUnitSquare::LeafIntersectionIterator IntersectionIterator;
//             typedef ALUCubeUnitSquare::LeafGridView GV;
//             typedef ALUCubeUnitSquare::ctype ctype;

//             // get view
//             const GV& gv=grid.leafView();

//             // Do some random refinement. The result is a grid that may
//             // contain multiple hanging nodes per edge.
//             for(int i=0; i<4;++i){
//                 Iterator it = grid.leafbegin<0,Dune::All_Partition>();
//                 Iterator eit = grid.leafend<0,Dune::All_Partition>();

//                 for(;it!=eit;++it){
//                     if((double)rand()/(double)RAND_MAX > 0.6)
//                         grid.mark(1,*(it));
//                 }
//                 grid.preAdapt();
//                 grid.adapt();
//                 grid.postAdapt();
//             }

//             // solve problem
//             mimetictest(gv,"mimetic_ALUCUBE_3d");
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
