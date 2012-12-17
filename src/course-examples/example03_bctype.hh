/** \brief constraint parameter class selecting boundary condition type */


class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters 
{
  
  double time;

public:

  template<typename I>
  bool isDirichlet(
                   const I & intersection,   
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    
    if( xg[0]>1.0-1E-6 )
      return false; // no Dirichlet b.c. on the eastern boundary
    
    return true;  // Dirichlet b.c. on all other boundaries
  }


  //! set time for subsequent evaluation
  void setTime (double t) 
  {
    time = t;
  }

};
