// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
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
#include<dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include<dune/pdelab/localoperator/diffusion.hh>

#include"gridexamples.hh"
#include"permeability_generator.hh"

//===============================================================
// Parameter functions for single phase flow problem
//===============================================================

// function for visualizing a random permeability field
template<typename GV, typename RF>
class RandomPermeability
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  RandomPermeability<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,RandomPermeability<GV,RF> > BaseT;

  RandomPermeability (const GV& gv, Dune::FieldVector<double,GV::dimension> correlation_length,
					  double variance = 1.0, double mean = 0.0, long modes = 1000, long seed = -1083) 
	: BaseT(gv), field(correlation_length,variance,mean,modes,seed)
  {
  }

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = log10(field.eval(x));
  }
private:
  EberhardPermeabilityGenerator<GV::dimension> field;
};


// function for defining the diffusion tensor
template<typename GV, typename RF>
class K_D
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_D<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_D<GV,RF> > BaseT;

  K_D (const GV& gv_, Dune::FieldVector<double,GV::dimension> correlation_length,
	   double variance = 1.0, double mean = 0.0, long modes = 1000, long seed = -1083) 
	: gv(gv_), is(gv.indexSet()), perm(is.size(0))
  {
	typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
	typedef typename Traits::DomainFieldType DF;
	const int dim = GV::dimension;
	double mink=1E100;
	double maxk=-1E100;

	
	EberhardPermeabilityGenerator<GV::dimension> field(correlation_length,variance,mean,modes,seed);

	for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
	  {
		int id = is.index(*it);
        Dune::GeometryType gt = it->geometry().type();
        Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
        Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
		perm[id]=field.eval(globalcenter);
		mink = std::min(mink,log10(perm[id]));
		maxk = std::max(maxk,log10(perm[id]));
	  }
	std::cout << "log10(mink)=" << mink << " log10(maxk)=" << maxk << std::endl;
  }

  K_D ( const GV& gv_, const std::vector<RF>& perm_)
    : gv(gv_), is(gv.indexSet()), perm(perm_)
  {}
    
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  { 
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          y[i][i] = perm[is.index(e)];
        else
          y[i][j] = 0.0;
  }
  
  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

  inline const RF& getElementPermeability(const typename GV::template Codim<0>::EntityPointer& e) const
  {
    return perm[is.index(*e)];
  }
  
private:
  const GV& gv;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;
};

// function for defining the source term
template<typename GV, typename RF>
class A0_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  A0_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,A0_D<GV,RF> > BaseT;

  A0_D (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class F_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_D<GV,RF> > BaseT;

  F_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0; 
  }
};

// boundary grid function selecting boundary conditions 
template<typename GV>
class B_D
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<GV,int,1,
                                                                             Dune::FieldVector<int,1> >,
                                                  B_D<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_D<GV> > BaseT;

  B_D (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig, 
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {  
	Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      xg = ig.geometry().global(x);

    if (xg[0]<1E-6 || xg[0]>1.0-1E-6)
	  y = 1; // Dirichlet
	else
	  y = 0;
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G_D<GV,RF> > BaseT;

  G_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0;
    if (x[0]<1E-6)
	  y = 1;
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J_D
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J_D<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J_D<GV,RF> > BaseT;

  J_D (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x, 
							  typename Traits::RangeType& y) const
  {
	y = 0;
	return;
  }
};

//===============================================================
// set up diffusion problem and solve it
//===============================================================

template<typename BType, typename GType, typename KType, typename A0Type, typename FType, typename JType,
         typename GV, typename FEM> 
void driver (const BType& b, const GType& g, 
             const KType& k, const A0Type& a0, const FType& f, const JType& j,
             const GV& gv, const FEM& fem, std::string filename)
{
  // constants and types and global variables
  typedef typename GV::Grid::ctype DF;
  typedef typename FEM::Traits::LocalFiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;
  Dune::Timer watch;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,Dune::PDELab::ConformingDirichletConstraints,
    Dune::PDELab::ISTLVectorBackend<1> > GFS; 
  GFS gfs(gv,fem);

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(b,gfs,cg);

  // make coefficent Vector and initialize it from a function
  typedef typename GFS::template VectorContainer<R>::Type V;
  V x(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,x);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x);

  // make grid function operator
  typedef Dune::PDELab::Diffusion<KType,A0Type,FType,BType,JType,2> LOP; 
  LOP lop(k,a0,f,b,j);
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LOP,C,C,Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,cg,gfs,cg,lop);

  // represent operator as a matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  watch.reset();
  M m(gos); m = 0.0;
  std::cout << "=== matrix setup " <<  watch.elapsed() << " s" << std::endl;
  watch.reset();
  gos.jacobian(x,m);
  std::cout << "=== jacobian assembly " <<  watch.elapsed() << " s" << std::endl;
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // solve the jacobian system
  V r(gfs,0.0);
  gos.residual(x,r);
  V z(gfs,0.0);

  // set up parallel solver
  typedef Dune::PDELab::NonoverlappingHelper<GFS> PHELPER;
  PHELPER phelper(gfs);
  typedef Dune::PDELab::NonoverlappingOperator<GFS,M,V,V> POP;
  POP pop(gfs,m,phelper);
  typedef Dune::PDELab::NonoverlappingScalarProduct<GFS,V> PSP;
  PSP psp(gfs,phelper);
  typedef Dune::PDELab::NonoverlappingRichardson<GFS,V,V> PRICH;
  PRICH prich(gfs,phelper);
  int verbose;
  if (gv.comm().rank()==0) verbose=1; else verbose=0;
  Dune::CGSolver<V> solver(pop,psp,prich,1E-8,40000,verbose);
  Dune::InverseOperatorResult stat;  
  solver.apply(z,r,stat);
  x -= z;

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunction<GFS,V> DGF;
  DGF dgf(gfs,x);
  
  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(dgf,"solution"));
  vtkwriter.pwrite(filename.c_str(),"vtk","",Dune::VTKOptions::binaryappended);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
		if(helper.rank()==0)
		  std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    
#if HAVE_MPI
    // Q1, 2d
    if (false)
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(4);
      Dune::FieldVector<bool,2> B(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,B,overlap);
      grid.globalRefine(4);
      
      typedef Dune::YaspGrid<2>::ctype DF;
      typedef Dune::PDELab::Q12DLocalFiniteElementMap<DF,double> FEM;
      FEM fem;

      Dune::FieldVector<double,2> correlation_length;
      correlation_length = 1.0/64.0;
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafView();
      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_yasp2d_Q1");
    }
#endif

#if HAVE_MPI
    // Q1, 3d
    if (true)
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::FieldVector<int,3> N(4);
      Dune::FieldVector<bool,3> B(false);
      int overlap=0;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,B,overlap);
      grid.globalRefine(4);

      typedef Dune::YaspGrid<3>::ctype DF;
      typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,double,3> FEM;
      FEM fem;

      Dune::FieldVector<double,3> correlation_length;
      correlation_length = 1.0/64.0;
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafView();
      driver(B_D<GV>(gv), G_D<GV,double>(gv),
             K_D<GV,double>(gv,correlation_length,1.0,0.0,5000,-1083),
             A0_D<GV,double>(gv),F_D<GV,double>(gv),J_D<GV,double>(gv),
             gv,fem,"single_phase_yasp3d_Q1");
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
