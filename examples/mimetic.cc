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
#include<dune/pdelab/finiteelementmap/p12dfem.hh>
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>
#include<dune/pdelab/finiteelementmap/pk3dfem.hh>
#include<dune/pdelab/finiteelementmap/q12dfem.hh>
#include<dune/pdelab/finiteelementmap/q22dfem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/finiteelementmap/hangingnodeconstraints.hh>
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
#include<dune/pdelab/localoperator/poisson.hh>

#include "gridexamples.hh"

/* ==============================================

   List of changes to other modules

   finiteelement/common/localcoefficents.hh -> intersectionCodim tag
   dune/pdelab/gridfunctionspace/gridfunctionspace.hh -> new gridfunctionspace


   ==============================================
 */

//===============================================================
// Index set for intersections
// implementation only for 
//  - hanging node grids, but conforming macro grid
//  - no multiple element types !
//===============================================================

#include<vector>
#include<dune/common/exceptions.hh>

template<typename GV>
class IntersectionIndexSet
{
public:
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;
  typedef typename GV::Intersection Intersection;
  typedef typename GV::Traits::template Codim<0>::Entity Element; 
 
  IntersectionIndexSet (const GV& gv_)
    : gv(gv_), is(gv.indexSet())
  {
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator; 

    pre();
	for (ElementIterator it = gv.template begin<0>(); 
		 it!=gv.template end<0>(); ++it)
      {
        visit(*it);
      }
    post();
  }

  // number of intersections in index set 
  // (intersections do not have a geometry type)
  IndexType size () const
  {
    return intersection_counter;
  }

  // number of intersections associated with given element
  IndexType size (const Element& element) const
  {
    return number_of_intersections[is.index(element)];
  }

  // get index assigned to intersection
  IndexType index (const Intersection& intersection) const
  {
    // on the boundary, take own codim 1 subentity
    if (intersection.boundary() && (!intersection.neighbor()))
      return codim1index_to_intersectionindex[is.subIndex(*(intersection.inside()),intersection.indexInInside(),1)];
    
    // if we have a neighbor, take higher level
    if (intersection.neighbor())
      {
        if (intersection.inside()->level()>=intersection.outside()->level())
          return codim1index_to_intersectionindex[is.subIndex(*(intersection.inside()),intersection.indexInInside(),1)];
        else
          return codim1index_to_intersectionindex[is.subIndex(*(intersection.outside()),intersection.indexInOutside(),1)];
      }

    // we are at a processor boundary
    DUNE_THROW(Dune::Exception,"intersection index at processor boundary requested");
  }

  // get index of i'th intersection of element 
  // (in order they are visited by intersection iterator)
  IndexType subIndex (const Element& element, int i)
  {
    return element_intersection_subindex[entry[is.index(element)]+i];
  }

private:

  // prepare loop over elements
  void pre ()
  {
    codim1index_to_intersectionindex.resize(gv.size(1));
    invalidIndex = gv.size(1)+1;
    for (size_t i=0; i<codim1index_to_intersectionindex.size(); ++i)
      codim1index_to_intersectionindex[i] = invalidIndex;
    number_of_intersections.resize(gv.size(0));
    entry.resize(gv.size(0));
    element_intersection_subindex.resize(2*gv.size(1));
    intersection_counter = 0;
    oriented_intersection_counter = 0;
    std::cout << "number of codim 1 entities is " << gv.size(1) << std::endl;
  }

  // process given element
  void visit (const Element& element)
  {
    typedef typename GV::IntersectionIterator IntersectionIterator;
    size_t count = 0;
    entry[is.index(element)] = oriented_intersection_counter;
    IntersectionIterator endit = gv.iend(element);
    for (IntersectionIterator iit = gv.ibegin(element); iit!=endit; ++iit)
      {
        if (iit->neighbor())
          {
            IndexType c1index;
            if (iit->inside()->level()>=iit->outside()->level())
              c1index = is.subIndex(*(iit->inside()),iit->indexInInside(),1);
            else
              c1index = is.subIndex(*(iit->outside()),iit->indexInOutside(),1);
            if (codim1index_to_intersectionindex[c1index]==invalidIndex)
              codim1index_to_intersectionindex[c1index]=intersection_counter++;
            element_intersection_subindex[oriented_intersection_counter] = codim1index_to_intersectionindex[c1index];
          }
        else if (iit->boundary())
          {
            IndexType c1index = is.subIndex(*(iit->inside()),iit->indexInInside(),1);
            if (codim1index_to_intersectionindex[c1index]==invalidIndex)
              codim1index_to_intersectionindex[c1index]=intersection_counter++;
            element_intersection_subindex[oriented_intersection_counter] = codim1index_to_intersectionindex[c1index];
         }
        count++;
        oriented_intersection_counter++;
      }
    number_of_intersections[is.index(element)] = static_cast<unsigned char>(count);
  }

  // finalize computation of index set
  void post ()
  {
    std::cout << "number of oriented intersections " << oriented_intersection_counter << std::endl;
    std::cout << "number of intersections " << intersection_counter << std::endl;
  }

  const GV& gv;       // our grid view
  const IndexSet& is; // index set of the grid view

  // we assume that the mesh is conforming in the
  // sense that there is a codim 1 entity for each intersection
  // the following vector assigns intersection to the codim 1 entity on the "higher" side
  std::vector<IndexType> codim1index_to_intersectionindex;

  // number of intersections of an element
  std::vector<unsigned char> number_of_intersections;

  // map (element index,interection number) to index 
  std::vector<IndexType> element_intersection_subindex;

  // entry point of element into element_intersection_subindex
  std::vector<size_t> entry;

  size_t intersection_counter;
  size_t oriented_intersection_counter;
  IndexType invalidIndex;
};

//===============================================================
// implement a local finite element and map for mimetic ...
//===============================================================

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/finiteelements/common/localbasis.hh>
#include<dune/finiteelements/common/localinterpolation.hh>
#include<dune/finiteelements/common/localcoefficients.hh>
#include<dune/finiteelements/common/localfiniteelement.hh>
#include<dune/pdelab/finiteelementmap/finiteelementmap.hh>

template<class D, class R, int dim>
class MimeticLocalBasis : 
  public Dune::C0LocalBasisInterface<
  Dune::C0LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,
						   R,1,Dune::FieldVector<R,1> >,
  MimeticLocalBasis<D,R,dim> > 
{
public:
  typedef Dune::C0LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,
                                   R,1,Dune::FieldVector<R,1> > Traits;

  MimeticLocalBasis (unsigned int variant_)
    : variant(variant_)
  {
  }

  MimeticLocalBasis ()
    : variant(0)
  {
  }

  unsigned int size () const { return variant; }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (
				 const typename Traits::DomainType& in,
		   std::vector<typename Traits::RangeType>& out) const 
  { 
    DUNE_THROW(Dune::Exception,"mimetic basis evaluation not available");    
  }
  
  //! \brief Polynomial order of the shape functions
  unsigned int order () const  
  {
    DUNE_THROW(Dune::Exception,"mimetic order evaluation not available");    
  }

private:
  unsigned int variant;
};

template<class LB>
class MimeticLocalInterpolation : public Dune::
  LocalInterpolationInterface<MimeticLocalInterpolation<LB> > {
public:

  //! \brief Local interpolation of a function
  template<typename F, typename C>
  void interpolate (const F& f, std::vector<C>& out) const {
    DUNE_THROW(Dune::Exception,"mimetic local interpolation not available");    
  }
};

class MimeticLocalCoefficients 
  : public Dune::LocalCoefficientsInterface<MimeticLocalCoefficients> {
public:
  MimeticLocalCoefficients (unsigned int variant_) 
    : variant(variant_), li(variant_)  
  {
	for (int i=0; i<variant; i++) li[i] = Dune::LocalKey(i,Dune::intersectionCodim,0);
  }
  
  MimeticLocalCoefficients () 
    : variant(0), li(0)  
  {
  }
  
  //! number of coefficients
  int size () const { return variant; }
  
  //! map index i to local key
  const Dune::LocalKey& localKey (int i) const  {
	return li[i];
  }
  
private:
  unsigned int variant;
  std::vector<Dune::LocalKey> li;
};

template<class D, class R, int dim>
class MimeticLocalFiniteElement 
  : public Dune::LocalFiniteElementInterface< 
  Dune::LocalFiniteElementTraits<MimeticLocalBasis<D,R,dim>,
                                 MimeticLocalCoefficients,
					  MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > >, 
  MimeticLocalFiniteElement<D,R,dim> >
{
  Dune::GeometryType gt;
  MimeticLocalBasis<D,R,dim> basis;
  MimeticLocalCoefficients coefficients;
  MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > interpolation;

public:
  typedef Dune::LocalFiniteElementTraits<MimeticLocalBasis<D,R,dim>,
                                         MimeticLocalCoefficients,
                MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > > Traits;

  MimeticLocalFiniteElement () 
  {}

  MimeticLocalFiniteElement (Dune::GeometryType::BasicType basicType) 
    : gt(basicType,dim)
  {}

  MimeticLocalFiniteElement (Dune::GeometryType::BasicType basicType, unsigned int variant) 
    : gt(basicType,dim), basis(variant), coefficients(variant)
  {}

  const typename Traits::LocalBasisType& localBasis () const 
  {
	return basis;
  }
  
  const typename Traits::LocalCoefficientsType& localCoefficients () const 
  {
	return coefficients;
  }
  
  const typename Traits::LocalInterpolationType& localInterpolation () const 
  {
	return interpolation;
  }
	
  Dune::GeometryType type () const { return gt; }
};


template<typename IIS, typename D, typename R, int dim>
class MimeticLocalFiniteElementMap : 
  public Dune::PDELab::LocalFiniteElementMapInterface<Dune::PDELab::LocalFiniteElementMapTraits< MimeticLocalFiniteElement<D,R,dim> >, 
                                                      MimeticLocalFiniteElementMap<IIS,D,R,dim> >,
  public Dune::PDELab::Countable
{
  typedef MimeticLocalFiniteElement<D,R,dim> FE;
     
public:
  //! \brief export type of the signature
  typedef Dune::PDELab::LocalFiniteElementMapTraits<FE> Traits;  

  //! \brief Use when Imp has a standard constructor
  MimeticLocalFiniteElementMap (const IIS& iis_, Dune::GeometryType::BasicType basicType) 
    : iis(iis_), bt(basicType)

  {
    // create a standard number of variants
    variant.resize(20);
    for (int i=0; i<20; i++) variant[i] = FE(bt,i);
  }

  //! \brief get local basis functions for entity
  template<class EntityType>
  const typename Traits::LocalFiniteElementType& find (const EntityType& e) const
  {
    size_t n = static_cast<size_t>(iis.size(e));
    if (n<variant.size())
      return variant[n];
    else
      {
        size_t old_n = variant.size();
        variant.resize(n+1);
        for (size_t i=old_n; i<n+1; i++) variant[i] = FE(bt,i);
        return variant[n];      
      }
  }

private:
  const IIS& iis;
  Dune::GeometryType::BasicType bt;
  mutable std::vector<FE> variant;
};


//===============================================================
// Problem setup and solution 
//===============================================================

// check DOFs in intersections
template<typename GV> 
void mimetictest (const GV& gv, std::string filename)
{
  // set up index set for intersections
  typedef IntersectionIndexSet<GV> IIS;
  IIS iis(gv);

  // make finite element map
  typedef MimeticLocalFiniteElementMap<IIS,double,double,GV::dimension> FEM;
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
