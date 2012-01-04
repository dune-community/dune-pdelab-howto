// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>
template<class U, class GFS, class X> 
double l2interpolationerror (const U& u, const GFS& gfs, X& x, 
							 int qorder=1) {
  // constants and types
  typedef typename GFS::Traits::GridViewType GV;
  const int dim = GV::Grid::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator 
	ElementIterator;
  typedef typename GFS::Traits::FiniteElementType::
	Traits::LocalBasisType::Traits FTraits;
  typedef typename FTraits::DomainFieldType D;
  typedef typename FTraits::RangeFieldType R;
  typedef typename FTraits::RangeType RangeType;
  
  // make local function space
  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  LFS lfs(gfs);                            /*@\label{l2int:lfs}@*/
  std::vector<R> xl(lfs.maxSize());        // local coefficients
  std::vector<RangeType> b(lfs.maxSize()); // shape function values
  
  // loop over grid view
  double sum = 0.0;
  for (ElementIterator eit = gfs.gridview().template begin<0>();
	   eit!=gfs.gridview().template end<0>(); ++eit)
	{
	  lfs.bind(*eit);      // bind local function space to element /*@\label{l2int:bind}@*/
	  lfs.vread(x,xl);     // get local coefficients
	  Dune::GeometryType gt = eit->geometry().type();
	  const Dune::QuadratureRule<D,dim>&  /*@\label{l2int:quad}@*/
		rule = Dune::QuadratureRules<D,dim>::rule(gt,qorder);

	  for (typename Dune::QuadratureRule<D,dim>::const_iterator 
			 qit=rule.begin(); 
		   qit!=rule.end(); ++qit)
		{
		  // evaluate finite element function at integration point
		  RangeType u_fe(0.0);
		  lfs.finiteElement().localBasis().evaluateFunction(
            qit->position(),b);
		  for (unsigned int i=0; i<lfs.size(); i++)
			u_fe.axpy(xl[i],b[i]);

		  // evaluate the given grid function at integration point
		  RangeType u_given;
		  u.evaluate(*eit,qit->position(),u_given);

		  // accumulate error
		  u_fe -= u_given;
		  sum += u_fe.two_norm2()*qit->weight()*
			eit->geometry().integrationElement(qit->position());
		}
	}
  return sqrt(sum);
}
