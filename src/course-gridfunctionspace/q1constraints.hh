class Q1Constraints {
public:
  enum{doBoundary=true};enum{doProcessor=false};
  enum{doSkeleton=false};enum{doVolume=false};

  template<typename B, typename I, typename LFS, typename T>
  void boundary (const B& b, const I& ig, const LFS& lfs, T& trafo) const
  {
    Dune::FieldVector<typename I::ctype,I::dimension-1>
      ip(0.5);                     // test edge midpoint
	bool isDirichlet =
      b.isDirichlet(ig,ip);        // eval condition type

	if (!isDirichlet) return;      // done

	typename T::RowType empty;     // need not interpolate
	if (ig.indexInInside()==0) { trafo[0]=empty; trafo[2]=empty; }
	if (ig.indexInInside()==1) { trafo[1]=empty; trafo[3]=empty; }
	if (ig.indexInInside()==2) { trafo[0]=empty; trafo[1]=empty; }
	if (ig.indexInInside()==3) { trafo[2]=empty; trafo[3]=empty; }
  }
};
