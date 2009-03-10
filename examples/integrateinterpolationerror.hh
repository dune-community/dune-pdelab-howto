#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
template<class U, class GFS, class X> 
double integrateinterpolationerror (const U& u, const GFS& gfs, X& x, int qorder=1)
{
  // constants and types
  typedef typename GFS::Traits::GridViewType GV;
  const int dim = GV::Grid::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GFS::Traits::LocalFiniteElementType::
	Traits::LocalBasisType::Traits FTraits;
  typedef typename FTraits::DomainFieldType D;
  typedef typename FTraits::RangeFieldType R;
  typedef typename FTraits::RangeType RangeType;
  
  // make local function space
  typedef typename GFS::LocalFunctionSpace LFS;
  LFS lfs(gfs);
  std::vector<R> xl(lfs.maxSize());        // stores local coefficients
  std::vector<RangeType> b(lfs.maxSize()); // stores shape function values
  
  // loop over grid view
  double sum = 0.0;
  for (ElementIterator eit = gfs.gridview().template begin<0>();
	   eit!=gfs.gridview().template end<0>(); ++eit)
	{
	  lfs.bind(*eit);      // bind local function space to element
	  lfs.vread(x,xl);     // get local coefficients
	  Dune::GeometryType gt = eit->geometry().type();
	  const Dune::QuadratureRule<D,dim>& 
		rule = Dune::QuadratureRules<D,dim>::rule(gt,qorder);

	  for (typename Dune::QuadratureRule<D,dim>::const_iterator qit=rule.begin(); 
		   qit!=rule.end(); ++qit)
		{
		  // evaluate finite element function at integration point
		  RangeType u_fe(0.0);
		  lfs.localFiniteElement().localBasis().evaluateFunction(qit->position(),b);
		  for (int i=0; i<lfs.size(); i++)
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
