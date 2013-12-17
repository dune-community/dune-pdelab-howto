// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Demonstrate parallel data exchange with a gridfunctionspace 
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

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/tags.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/localoperator/laplacedirichletccfv.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

//==============================================================================
// driver
//==============================================================================
int rank;

class TestGatherScatter
{
public:
  TestGatherScatter ()
  {
    gather_n = scatter_n = 0;
  }

  template<class MessageBuffer, class DataType>
  void gather (MessageBuffer& buff, DataType& data)
  {
    buff.write(data);
    std::cout << "[" << rank << "] gather: count=" << gather_n << " adress=" << &data << " data=" << data << std::endl;
    gather_n++;
  }
  
  template<class MessageBuffer, class DataType>
  void scatter (MessageBuffer& buff, DataType& data)
  {
    DataType x; 
    buff.read(x);
    std::cout << "[" << rank << "] scatter: count=" << scatter_n << " received=" << x << " adress=" << &data << " data=" << data << std::endl;
    data = x;
    scatter_n++;
  }
private:
  int gather_n;
  int scatter_n;
};


template<class GV> 
void testp0 (const GV& gv)
{
  // some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
  
  // make grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::P0ParallelConstraints,
    Dune::PDELab::ISTLVectorBackend<>
    > GFS;
  watch.reset();
  Dune::PDELab::P0ParallelConstraints con;
  GFS gfs(gv,fem,con);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // make vector for old time step and initialize
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
  V v(gfs);

  // test generic communication
  typename V::iterator vi = v.begin();
  for (size_t i=0; vi!=v.end(); i++, ++vi) *vi = 0.1*i;
    
  // test default data handles
  Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,v);
  gv.communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::CopyDataHandle<GFS,V> copydh(gfs,v);
  gv.communicate(copydh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::MinDataHandle<GFS,V> mindh(gfs,v);
  gv.communicate(mindh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::MaxDataHandle<GFS,V> maxdh(gfs,v);
  gv.communicate(maxdh,Dune::All_All_Interface,Dune::ForwardCommunication);

  // test power
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,3,
    Dune::PDELab::ISTLVectorBackend<>,
    Dune::PDELab::LexicographicOrderingTag
  > PGFS;
  PGFS pgfs(gfs);
  typedef typename Dune::PDELab::BackendVectorSelector<PGFS,RF>::Type PV;
  PV pv(pgfs);
  typename PV::iterator pvi = pv.begin();
  for (size_t i=0; pvi!=pv.end(); i++, ++pvi) *pvi = i;

  // test composite
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    Dune::PDELab::ISTLVectorBackend<>,
      Dune::PDELab::LexicographicOrderingTag,
	  GFS,PGFS> CGFS;              
  CGFS cgfs(gfs,pgfs);
  
  typedef typename Dune::PDELab::BackendVectorSelector<CGFS,RF>::Type CV;
  CV cv(cgfs);
  typename PV::iterator cvi = cv.begin();
  for (size_t i=0; cvi!=cv.end(); i++, ++cvi) *cvi = 100*i;
  Dune::PDELab::AddDataHandle<CGFS,CV> adddh2(cgfs,cv);
  gv.communicate(adddh2,Dune::All_All_Interface,Dune::ForwardCommunication);
}

template<class GV> 
void testp1 (const GV& gv)
{
  // some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  Dune::Timer watch;

  // instantiate finite element maps
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,1> FEM;
  FEM fem(gv);
  
  // make grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,
    Dune::PDELab::ConformingDirichletConstraints,Dune::PDELab::ISTLVectorBackend<>
    > GFS;
  watch.reset();
  Dune::PDELab::ConformingDirichletConstraints con;
  GFS gfs(gv,fem,con);
  std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

  // make vector for old time step and initialize
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
  V v(gfs);

  // test generic communication
  typename V::iterator vi = v.begin();
  for (size_t i=0; vi!=v.end(); i++, ++vi) *vi = 0.1*i;

  // test default data handles
  Dune::PDELab::AddDataHandle<GFS,V> adh(gfs,v);
  gv.communicate(adh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::CopyDataHandle<GFS,V> cdh(gfs,v);
  gv.communicate(cdh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::MinDataHandle<GFS,V> mindh(gfs,v);
  gv.communicate(mindh,Dune::All_All_Interface,Dune::ForwardCommunication);
  Dune::PDELab::MaxDataHandle<GFS,V> maxdh(gfs,v);
  gv.communicate(maxdh,Dune::All_All_Interface,Dune::ForwardCommunication);

  // test power
  typedef Dune::PDELab::ISTLVectorBackend
    <Dune::PDELab::ISTLParameters::static_blocking,3> VBE;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,3,VBE,
    Dune::PDELab::EntityBlockedOrderingTag
    > PGFS;
  PGFS pgfs(gfs);
  typedef typename Dune::PDELab::BackendVectorSelector<PGFS,RF>::Type PV;
  PV pv(pgfs);
  typename PV::iterator pvi = pv.begin();
  for (size_t i=0; pvi!=pv.end(); i++, ++pvi) *pvi = i;
  Dune::PDELab::AddDataHandle<PGFS,PV> adh2(pgfs,pv);
  gv.communicate(adh2,Dune::All_All_Interface,Dune::ForwardCommunication);
}

//==============================================================================
// grid setup
//==============================================================================

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
		  std::cout << "parallel run on " << helper.size() << " processes" << std::endl;
	  }
    rank = helper.rank();

    // 2D
    if (true)
    {
      // make grid
      Dune::FieldVector<double,2> L; L[0] = 2.0; L[1] = 1.0;
      Dune::array<int,2> N;
      N[0] = 4;   N[1] = 2;
      std::bitset<2> periodic(false);
      int overlap=1;
      Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
 
      // solve problem :
      testp0(grid.leafGridView());
      testp1(grid.leafGridView());
    }

    // 3D
    if (false)
    {
      // make grid
      Dune::FieldVector<double,3> L; L[0] = 2; L[1] = 1; L[2] = 1;
      Dune::array<int,3> N;
      N[0] = 4;    N[1] = 2;    N[2] = 2;
      std::bitset<3> periodic(false);
      int overlap=1;
      Dune::YaspGrid<3> grid(helper.getCommunicator(),L,N,periodic,overlap);
      
      // solve problem :)
      testp0(grid.leafGridView());
      testp1(grid.leafGridView());
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
