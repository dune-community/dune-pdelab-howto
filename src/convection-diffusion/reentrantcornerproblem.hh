// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_REENTRANTCORNERPROBLEM_HH
#define DUNE_PDELAB_REENTRANTCORNERPROBLEM_HH

template<typename GV, typename RF>
class ReentrantCornerProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1.0 : 0.0;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType 
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType 
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    //typename Traits::DomainType xglobal = is.geometry().global(x);
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType 
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);

    if( x[0]>-1e-12 && x[1]<1e-12 )
        return 0;
    else{
      typename Traits::DomainFieldType theta = std::atan2(x[1], x[0]);
      if(theta < 0.0) theta += 2*M_PI;
      typename Traits::DomainFieldType r = x.two_norm();
      return pow(r,2.0/3.0)*std::sin(theta*2.0/3.0);
    }
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType 
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }
};

//! exact gradient of solution
template<typename GV, typename RF>
class ExactGradient
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
                                                                           GV::dimension,Dune::FieldVector<RF,GV::dimension> >,
                                          ExactGradient<GV,RF> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
                                           GV::dimension,Dune::FieldVector<RF,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,ExactGradient<GV,RF> > BaseT;

  ExactGradient (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    typename Traits::DomainType x = e.geometry().global(xlocal);

    typename Traits::DomainFieldType theta = std::atan2(x[1], x[0]);
    if(theta < 0.0) theta += 2*M_PI;
    typename Traits::DomainFieldType r = x.two_norm();

    y[0] = (2.0/3.0)*pow(r,-1.0/3.0)*(cos(theta)*sin(2.0*theta/3.0) - sin(theta)*cos(2.0*theta/3.0));
    y[1] = (2.0/3.0)*pow(r,-1.0/3.0)*(sin(theta)*sin(2.0*theta/3.0) + cos(theta)*cos(2.0*theta/3.0));
  }
  
  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
};

/*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceAdapter 
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,
                                   Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >, 
  DifferenceAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor 
  DifferenceAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x, 
                        typename Traits::RangeType& y) const
  {  
    typename Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }
  
private:
  const T1& t1;
  const T2& t2;
};

/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter 
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,
                                   Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >,
  DifferenceSquaredAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor 
  DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& x, 
                        typename Traits::RangeType& y) const
  {  
    typename T1::Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename T2::Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1.two_norm2();
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }
  
private:
  const T1& t1;
  const T2& t2;
};


#endif
