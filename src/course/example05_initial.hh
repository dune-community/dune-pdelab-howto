#include"../utility/permeability_generator.hh"

/** \brief A function for initial values of u_0
 */
template<typename GV, typename RF>
class U0Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U0Initial<GV,RF> >
{
  const GV& gv;
  Dune::FieldVector<double,GV::dimension> correlation_length;
  EberhardPermeabilityGenerator<GV::dimension> field;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  U0Initial (const GV& gv_) 
    : gv(gv_), correlation_length(2.0/100.0), field(correlation_length,1,0.0,5000,-1083) 
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    y = log10(field.eval(x));
    return;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of u_1
 */
template<typename GV, typename RF>
class U1Initial
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, U1Initial<GV,RF> >
{
  const GV& gv;
  Dune::FieldVector<double,GV::dimension> correlation_length;
  EberhardPermeabilityGenerator<GV::dimension> field;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  U1Initial (const GV& gv_) 
    : gv(gv_), correlation_length(2.0/100.0), field(correlation_length,1,0.0,5000,-34) 
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    y = log10(field.eval(x));
    return;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
