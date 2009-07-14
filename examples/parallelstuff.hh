#ifndef PARALLEL_STUFF_HH
#define PARALLEL_STUFF_HH

#include<dune/istl/solvercategory.hh>
#include<dune/istl/scalarproducts.hh>
#include<dune/istl/operators.hh>

//==============================================================
// partition all entities of a given codim
// assumes they all have the same geometry type
//==============================================================
template<class GV, int CODIM>
class EntityPartitioning
{
  // A DataHandle class to build minimum of a std::vector
  // which is accompanied by an index set
  template<class IS, class V> // mapper type and vector type
  class MinimumExchange
	: public Dune::CommDataHandleIF<MinimumExchange<IS,V>,typename V::value_type>
  {
  public:
	//! export type of data for message buffer
	typedef typename V::value_type DataType;

	//! returns true if data for this codim should be communicated
	bool contains (int dim, int codim) const
	{
	  return (codim==CODIM);
	}

	//! returns true if size per entity of given dim and codim is a constant
	bool fixedsize (int dim, int codim) const
	{
	  return true;
	}

	/*! how many objects of type DataType have to be sent for a given entity

	  Note: Only the sender side needs to know this size. 
	*/
	template<class EntityType>
	size_t size (EntityType& e) const
	{
	  return 1;
	}

	//! pack data from user to message buffer
	template<class MessageBuffer, class EntityType>
	void gather (MessageBuffer& buff, const EntityType& e) const
	{
	  buff.write(v[indexset.index(e)]);
	}

	/*! unpack data from message buffer to user

	  n is the number of objects sent by the sender
	*/
	template<class MessageBuffer, class EntityType>
	void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	{
	  DataType x; buff.read(x);
	  if (x>=0) // other is -1 means, he does not want it 
		v[indexset.index(e)] = std::min(x,v[indexset.index(e)]);
	}

	//! constructor
	MinimumExchange (const IS& indexset_, V& v_) 
	  : indexset(indexset_), v(v_)
	{}
 
  private:
	const IS& indexset;
	V& v;
  };

public:
  typedef typename GV::Traits::template Codim<CODIM>::Iterator Iterator; 
  typedef typename GV::Traits::template Codim<CODIM>::Entity Entity; 

  /*! \brief Constructor needs to know the grid function space
   */
  EntityPartitioning (const GV& gridview_) :
	gridview(gridview_), assignment(gridview_.size(CODIM)), rank(gridview.comm().rank())
  {
	// extract types
	typedef typename GV::IndexSet IndexSet;

	// assign own rank to entities that I might have
	for (Iterator it = gridview.template begin<CODIM>(); 
		 it!=gridview.template end<CODIM>(); ++it)
	  if (it->partitionType()==Dune::InteriorEntity || it->partitionType()==Dune::BorderEntity)
		assignment[gridview.indexSet().template index(*it)] = rank; // set to own rank
	  else
		assignment[gridview.indexSet().template index(*it)] = -1; // I will not possibly own that vertex

	// exchange
	if (CODIM>0) // for CODIM==0 no communication is necessary
	  {
		MinimumExchange<IndexSet,std::vector<int> > dh(gridview.indexSet(),assignment);
		gridview.communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);
	  }

	// convert vector of minimum ranks to assignment vector
	for (Iterator it = gridview.template begin<CODIM>(); 
		 it!=gridview.template end<CODIM>(); ++it)
	  if (assignment[gridview.indexSet().template index(*it)]==rank)
		assignment[gridview.indexSet().template index(*it)] = 1;
	  else
		assignment[gridview.indexSet().template index(*it)] = 0;
  }

  bool owner (const Entity& entity) const
  {
	return assignment[gridview.indexSet().template index(entity)];
  }

  bool owner (size_t index) const
  {
	return assignment[index];
  }

  template<typename X>
  void keepOwner (X& x) const
  {
	for (size_t i=0; i<x.size(); ++i) x[i] *= assignment[i];
  }

private:
  const GV& gridview;
  std::vector<int> assignment;
  int rank;
};



//==============================================================
// scalar product for nonoverlapping decomposition
// works on dofs assigned to single codim
//==============================================================
template<class GV, int CODIM, class X>
class SimpleNonoverlappingScalarProduct : public Dune::ScalarProduct<X>
{
public:
  //! export types
  typedef X domain_type;
  typedef typename X::ElementType field_type;

  //! define the category
  enum {category=Dune::SolverCategory::nonoverlapping};

  /*! \brief Constructor needs to know the grid function space
   */
  SimpleNonoverlappingScalarProduct (const GV& gridview_) :
	gridview(gridview_), partitioning(gridview_)
  {
  }

  /*! \brief Dot product of two vectors. 
	It is assumed that the vectors are consistent on the interior+border
	partition.
  */
  virtual field_type dot (const X& x, const X& y)
  {
	// do local scalar product on unique partition
	field_type sum = 0;
	for (int i=0; i<x.size(); ++i)
	  sum += (x[i]*y[i])*partitioning.owner(i); // visit each index exactly once

	// do global communication
	return gridview.comm().sum(sum);
  }

  /*! \brief Norm of a right-hand side vector. 
	The vector must be consistent on the interior+border partition
  */
  virtual double norm (const X& x)
  {
	return sqrt(static_cast<double>(this->dot(x,x)));
  }

  void uniquefy (X& x) const
  {
	partitioning.keepOwner(x);
  }

private:
  const GV& gridview;
  EntityPartitioning<GV,CODIM> partitioning;
};


//==============================================================
// operator for nonoverlapping decomposition
// works on dofs assigned to single codim
//==============================================================
template<class GV, int CODIM, class M, class X, class Y>
class SimpleNonoverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
{
  // A DataHandle class to sum entries of P1 DOF vector
  // which is accompanied by an index set
  template<class IS, class V> // mapper type and vector type
  class AddingDataHandle
	: public Dune::CommDataHandleIF<AddingDataHandle<IS,V>,typename V::value_type>
  {
  public:
	//! export type of data for message buffer
	typedef typename V::block_type DataType;

	//! returns true if data for this codim should be communicated
	bool contains (int dim, int codim) const
	{
	  return (codim==CODIM);
	}

	//! returns true if size per entity of given dim and codim is a constant
	bool fixedsize (int dim, int codim) const
	{
	  return true;
	}

	/*! how many objects of type DataType have to be sent for a given entity

	  Note: Only the sender side needs to know this size. 
	*/
	template<class EntityType>
	size_t size (EntityType& e) const
	{
	  return 1;
	}

	//! pack data from user to message buffer
	template<class MessageBuffer, class EntityType>
	void gather (MessageBuffer& buff, const EntityType& e) const
	{
	  buff.write(v[indexset.index(e)]);
	}

	/*! unpack data from message buffer to user

	  n is the number of objects sent by the sender
	*/
	template<class MessageBuffer, class EntityType>
	void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	{
	  DataType x; buff.read(x);
	  v[indexset.index(e)] += x;
	}

	//! constructor
	AddingDataHandle (const IS& indexset_, V& v_) 
	  : indexset(indexset_), v(v_)
	{}
 
  private:
	const IS& indexset;
	V& v;
  };

public:
  //! export types
  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;

  //redefine the category, that is the only difference
  enum {category=Dune::SolverCategory::nonoverlapping};

  SimpleNonoverlappingOperator (const GV& gridview_, const M& A) 
	: gridview(gridview_), _A_(A), partitioning(gridview_) 
  {
  }

  //! apply operator to x:  \f$ y = A(x) \f$
  virtual void apply (const X& x, Y& y) const
  {
	// apply local operator; now we have sum y_p = sequential y
	_A_.mv(x,y);

	// accumulate y on border
	if (CODIM>0)
	  {
		AddingDataHandle<typename GV::IndexSet,Y> dh(gridview.indexSet(),y);
		gridview.communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }
	//	partitioning.keepOwner(y);
  }
  
  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
  {
	// apply local operator; now we have sum y_p = sequential y
	_A_.usmv(alpha,x,y);

	// accumulate y on border
	if (CODIM>0)
	  {
		AddingDataHandle<typename GV::IndexSet,Y> dh(gridview.indexSet(),y);
		gridview.communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }
	//	partitioning.keepOwner(y);
  }
  
  //! get matrix via *
  virtual const M& getmat () const
  {
	return _A_;
  }
  
private:
  const GV& gridview;
  const M& _A_;
  EntityPartitioning<GV,CODIM> partitioning;
};



//==============================================================
// data handle for communication in preconditioner
//==============================================================
template<class GV, int CODIM, class V> // mapper type and vector type
class AddingDataHandle
  : public Dune::CommDataHandleIF<AddingDataHandle<GV,CODIM,V>,typename V::ElementType>
{
  typedef typename GV::IndexSet IS;

public:
  //! export type of data for message buffer
  typedef typename V::ElementType DataType;

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
	return (codim==CODIM);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
	return true;
  }

  /*! how many objects of type DataType have to be sent for a given entity

	Note: Only the sender side needs to know this size. 
  */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
	return 2;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
	DataType x;
	x = v[indexset.index(e)][0]; 
	buff.write(x);
	x = v[indexset.index(e)][1]; 
	buff.write(x);
  }

  /*! unpack data from message buffer to user

	n is the number of objects sent by the sender
  */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
	DataType x; 
	buff.read(x);
	v[indexset.index(e)][0] += x;
	buff.read(x);
	v[indexset.index(e)][1] += x;
  }

  //! constructor
  AddingDataHandle (const GV& gv, V& v_) 
	: indexset(gv.indexSet()), v(v_)
  {}
 
private:
  const IS& indexset;
  V& v;
};



//==============================================================
// parallel Richardson preconditioner
//==============================================================
template<class GV, int CODIM, class X, class Y>
class NonoverlappingRichardson : public Dune::Preconditioner<X,Y> 
{
public:
  //! \brief The domain type of the preconditioner.
  typedef X domain_type;
  //! \brief The range type of the preconditioner.
  typedef Y range_type;
  //! \brief The field type of the preconditioner.
  typedef typename X::ElementType field_type;

  // define the category
  enum {
	//! \brief The category the preconditioner is part of.
	category=Dune::SolverCategory::nonoverlapping
  };

  /*! \brief Constructor.
      
    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    \param n The number of iterations to perform.
    \param w The relaxation factor.
  */
  NonoverlappingRichardson (const GV& gridview_, field_type w=1.0)
	: gridview(gridview_)
  {
	_w = w;
  }

  /*!
	\brief Prepare the preconditioner.
      
	\copydoc Preconditioner::pre(X&,Y&)
  */
  virtual void pre (X& x, Y& b) {}

  /*!
	\brief Apply the precondioner.
      
	\copydoc Preconditioner::apply(X&,const Y&)
  */
  virtual void apply (X& v, const Y& d)
  {
	v = d;
	v *= _w;
	AddingDataHandle<GV,CODIM,X> dh(gridview,v);
	gridview.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  /*!
	\brief Clean up.
      
	\copydoc Preconditioner::post(X&)
  */
  virtual void post (X& x) {}

private:
  //! \brief The relaxation factor to use.
  field_type _w;
  const GV& gridview;
};


//==============================================================
// Block jacobi with any subdomain preconditioner
//==============================================================
template<class GV, int CODIM, class P, class C>
class NonoverlappingBlockJacobi 
  : public Dune::Preconditioner<typename P::domain_type,typename P::range_type> 
{
public:
  //! \brief The domain type of the preconditioner.
  typedef typename P::domain_type domain_type;
  //! \brief The range type of the preconditioner.
  typedef typename P::range_type range_type;

  // define the category
  enum {
	//! \brief The category the preconditioner is part of.
	category=Dune::SolverCategory::nonoverlapping
  };

  /*! \brief Constructor.
      
    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    \param n The number of iterations to perform.
    \param w The relaxation factor.
  */
  NonoverlappingBlockJacobi (const GV& gridview_, P& prec_, const C& c_)
	: gridview(gridview_), partitioning(gridview_), prec(prec_), c(c_)
  {
  }

  /*!
	\brief Prepare the preconditioner.
      
	\copydoc Preconditioner::pre(domain_type&,range_type&)
  */
  virtual void pre (domain_type& x, range_type& b) 
  {
	prec.pre(x,b);
  }

  /*!
	\brief Apply the precondioner.
      
	\copydoc Preconditioner::apply(domain_type&,const range_type&)
  */
  virtual void apply (domain_type& v, const range_type& d)
  {
	char buf[64];
	sprintf(buf,"[%02d] ",gridview.comm().rank());
	//Dune::printvector(std::cout,v.base(),"v on entry in preconditioner",buf,20,9,1);
	//Dune::printvector(std::cout,d.base(),"d on entry in preconditioner",buf,20,9,1);

	//------------------------------------------------------------------------
	// exchange defect because of missing matrix entries
	//------------------------------------------------------------------------
	range_type dd(d);
// 	partitioning.keepOwner(dd); // so it works even for overlap in codim 0
//  	AddingDataHandle<GV,CODIM,range_type> dh2(gridview,dd);
// 	gridview.communicate(dh2,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	set_constrained_dofs(c,0.0,dd);
	// end additional exchange
	//------------------------------------------------------------------------

	prec.apply(v,dd);
	partitioning.keepOwner(v); // so it works even for overlap in codim 0
 	AddingDataHandle<GV,CODIM,domain_type> dh(gridview,v);
	gridview.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	//Dune::printvector(std::cout,v.base(),"v after communication in preconditioner",buf,20,9,1);
  }

  /*!
	\brief Clean up.
      
	\copydoc Preconditioner::post(domain_type&)
  */
  virtual void post (domain_type& x) 
  {
	prec.post(x);
  }

private:
  const GV& gridview;
  EntityPartitioning<GV,CODIM> partitioning;
  P& prec;
  const C& c;
};


//==============================================================
// example setup
//==============================================================
/*
  SimpleNonoverlappingOperator<GV,CODIM,M,V,V> op(gv,m);
  SimpleNonoverlappingScalarProduct<GV,CODIM,V> sp(gv);
  NonoverlappingRichardson<GV,CODIM,V,V> prec(gv);
  int verbose;
  if (gv.comm().rank()==0) verbose=2; else verbose=0;
  Dune::CGSolver<V> solver(op,sp,prec,1E-8,20000,verbose);
*/

#endif
