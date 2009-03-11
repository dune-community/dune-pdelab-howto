class Q1Constraints
{
public:
  enum { doBoundary = true };
  enum { doSkeleton = false };
  enum { doVolume = false };

  template<typename B, typename I, typename LFS, typename T>
  void boundary (const B& b, const I& ig, const LFS& lfs, T& trafo) const
  {
	// 2D here, get midpoint of edge
	typename B::Traits::DomainType ip(0.5);

	// determine type of boundary condition
	typename B::Traits::RangeType bctype;
	b.evaluate(ig,ip,bctype);

	// >0 means Dirichlet boundary
	if (bctype<=0) return;

	// the two end nodes of the edge are constrained
	typename T::RowType empty;
	if (ig.numberInSelf()==0) { trafo[0] = empty; trafo[2] = empty; }
	if (ig.numberInSelf()==1) { trafo[1] = empty; trafo[3] = empty; }
	if (ig.numberInSelf()==2) { trafo[0] = empty; trafo[1] = empty; }
	if (ig.numberInSelf()==3) { trafo[2] = empty; trafo[3] = empty; }
  }
};
