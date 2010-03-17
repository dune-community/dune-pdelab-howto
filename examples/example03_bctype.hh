/** \brief boundary grid function selecting boundary conditions 0=Neumann, 1=Dirichlet */
template<typename GV>
class BCType : public Dune::PDELab::BoundaryGridFunctionBase<
Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> >,BCType<GV> >
{
  const GV& gv; double time;
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view
  BCType (const GV& gv_) : gv(gv_) {}

  //! return bc type at point on intersection
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const {  
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> x = i.geometry().global(xlocal);
    if (x[0]>1.0-1e-6) y = 0; else  y = 1; return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

  //! set time for subsequent evaluation
  void setTime (double t) {time = t;}
};
