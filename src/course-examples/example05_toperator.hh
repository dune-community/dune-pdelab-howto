#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

/** \brief A local operator for the mass operator (L_2 integral) in the system */
class Example05TimeLocalOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume<Example05TimeLocalOperator>,
    public Dune::PDELab::NumericalJacobianVolume<Example05TimeLocalOperator>,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // constructor remembers parameters
  Example05TimeLocalOperator (double tau_, unsigned int intorder_=2)
    : tau(tau_), intorder(intorder_) {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (assume Galerkin scheme U=V)
    typedef typename LFSU::template Child<0>::Type LFSU0;
    const LFSU0& lfsu0 = lfsu.template child<0>();
    typedef typename LFSU::template Child<1>::Type LFSU1;
    const LFSU1& lfsu1 = lfsu.template child<1>();

    // domain and range field type (assume both components have same RF)
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<RangeType> phi0(lfsu0.size());
        lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);
        std::vector<RangeType> phi1(lfsu1.size());
        lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);

        // compute u_0, u_1 at integration point
        RF u_0=0.0;
        for (size_type i=0; i<lfsu0.size(); i++)
	  u_0 += x(lfsu0,i)*phi0[i];
        RF u_1=0.0;
        for (size_type i=0; i<lfsu1.size(); i++)
	  u_1 += x(lfsu1,i)*phi1[i];

        // integration
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu0.size(); i++)
          r.accumulate(lfsu0,i,u_0*phi0[i]*factor);
        for (size_type i=0; i<lfsu1.size(); i++)
          r.accumulate(lfsu1,i,tau*u_1*phi1[i]*factor);
      }
  }
private:
  double tau;
  unsigned int intorder;
};
