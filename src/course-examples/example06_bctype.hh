/** \brief constraints parameter class selecting boundary conditions 
 */
template<typename GV>
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters 
{

  const GV& gv;

public:
  
  typedef GV GridViewType;

  //! construct from grid view
  BCTypeParam( const GV& gv_ ) 
    : gv(gv_) 
  {
  }
  
  template<typename I>
  bool isDirichlet(
                   const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    //Dune::FieldVector<typename I::ctype, I::dimension>
    //  xg = intersection.geometry().global( coord );

    return false; // Neumann b.c. everywhere
  }


  //! get a reference to the grid view
  inline const GV& getGridView () 
  {
    return gv;
  }

};
