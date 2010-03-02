/** \brief boundary grid function selecting boundary conditions 
 * 0 means Neumann
 * 1 means Dirichlet
 */
template<typename GV>
class BCType : public Dune::PDELab::BoundaryGridFunctionBase<
       Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
       Dune::FieldVector<int,1> >,BCType<GV> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view
  BCType (const GV& gv_) : gv(gv_) {}

  //! return bc type at point on intersection
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    Dune::FieldVector<typename GV::Grid::ctype,GV::dimension> 
      x = i.geometry().global(xlocal);
    y = 1; // default is Dirichlet
    if (x[0]>1.0-1e-6)
      y = 0; // Neumann
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
