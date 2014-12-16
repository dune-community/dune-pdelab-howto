/** \brief Boundary conditions for Hagen-Poiseuille flow.

    \tparam GV The underlying grid view.

    * The domain is either the unitsquare or the unitcube, depending on the dimension.

*/

template<typename GV>
class B_A
    : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                    BoundaryGridFunctionTraits<
                                                      GV, Dune::PDELab::StokesBoundaryCondition::Type, 1,
                                                      Dune::FieldVector<
                                                        Dune::PDELab::StokesBoundaryCondition::Type,1> >,
                                                    B_A<GV> >
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,BC::Type,1,Dune::FieldVector<BC::Type,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_A<GV> > BaseT;

  B_A (const GV& ) {}

  template<typename I>
  inline void evaluate (const I& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = BC::VelocityDirichlet;
  }
};

/** \brief Dirichlet velocity function for Hagen-Poiseuille flow.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.

    * Parabolic inflow velocity in tube pointing in x-direction.
    * The inflow velocity is normalized to 1.

*/

template<typename GV, typename RF>
class V_A
    : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                                                    V_A<GV,RF> >
{
public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_A<GV,RF> > BaseT;

    V_A (const GV& gv) : BaseT(gv) {}
    inline void evaluateGlobal (const typename Traits::DomainType& x,
        typename Traits::RangeType& y) const
    {
        y = 0;
        y[0] = x[1]*(1.0-x[1]);
        y[0] *= 4.0;
    }
};

/** \brief Neumann boundary function for Hagen-Poiseuille flow.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.

*/

template<typename GV, typename RF>
class P_A
    : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                    P_A<GV,RF> >
{
public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_A<GV,RF> > BaseT;

    P_A (const GV& gv) : BaseT(gv) {}
    inline void evaluateGlobal (const typename Traits::DomainType& x,
        typename Traits::RangeType& y) const
    {
        y = 0;
    }
};

/** \brief Source term function.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.

*/

template<typename GV, typename RF>
class F_A
    : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                                                    F_A<GV,RF> >
{
public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_A<GV,RF> > BaseT;

    F_A (const GV& gv) : BaseT(gv) {}
    inline void evaluateGlobal (const typename Traits::DomainType& x,
        typename Traits::RangeType& y) const
    {
        y = 0.0;
    }
};
