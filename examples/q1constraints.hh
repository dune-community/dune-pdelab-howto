class Q1Constraints
{
public:
  enum{doBoundary=true};enum{doSkeleton=false};enum{doVolume=false};

  template<typename B, typename I, typename LFS, typename T>
  void boundary (const B& b, const I& ig, const LFS& lfs, T& trafo) const
  {
	typename B::Traits::DomainType ip(0.5); // test edge midpoint
	typename B::Traits::RangeType bctype;   // return value
	b.evaluate(ig,ip,bctype);               // eval condition type

	if (bctype<=0) return;                  // done

	typename T::RowType empty;              // need not interpolate
	if (ig.numberInSelf()==0) { trafo[0]=empty; trafo[2]=empty; }
	if (ig.numberInSelf()==1) { trafo[1]=empty; trafo[3]=empty; }
	if (ig.numberInSelf()==2) { trafo[0]=empty; trafo[1]=empty; }
	if (ig.numberInSelf()==3) { trafo[2]=empty; trafo[3]=empty; }
  }
};
