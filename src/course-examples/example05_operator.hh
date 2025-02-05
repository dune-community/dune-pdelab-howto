#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

/** A local operator for solving the two-component reaction-diffusion
 *  system of Fitzhugh-Nagumo type. See also
 *  http://en.wikipedia.org/wiki/Reaction-diffusion
 *
 *   - d_0 \Delta u_0 - (\lambda*u_0 - u_0^3 - \sigma* u_1 + \kappa) = 0   in \Omega
 *   - d_1 \Delta u_1 - (u_0 - u_1)                                  = 0   in \Omega
 *
 * with natural boundary conditions
 *   \nabla u_0 \cdot n = 0   on \partial\Omega
 *   \nabla u_1 \cdot v = 0   on \partial\Omega
 *
 * with conforming finite elements on all types of grids in any dimension
 */
class Example05LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<Example05LocalOperator >,
  public Dune::PDELab::NumericalJacobianVolume<Example05LocalOperator >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // constructor stores parameters
  Example05LocalOperator (double d_0_, double d_1_, double lambda_, double sigma_,
                          double kappa_, unsigned int intorder_=2)
    : intorder(intorder_), d_0(d_0_), d_1(d_1_), lambda(lambda_),
      sigma(sigma_), kappa(kappa_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (assume Galerkin scheme U=V)
    typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
    const LFSU0& lfsu0 = lfsu.template child<0>();           // with template magic
    typedef typename LFSU::template Child<1>::Type LFSU1;
    const LFSU1& lfsu1 = lfsu.template child<1>();

    // domain and range field type (assume both components have same RF)
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator
           it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<RangeType> phi0(lfsu0.size());
        lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);
        std::vector<RangeType> phi1(lfsu1.size());
        lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);

        // compute u_0, u_1 at integration point
        RF u_0=0.0;
        for (size_type i=0; i<lfsu0.size(); i++)
          u_0 += x(lfsu0,i)*phi0[i];                // localIndex() maps dof within
        RF u_1=0.0;                                             // leaf space to all dofs
        for (size_type i=0; i<lfsu1.size(); i++)                // within given element
          u_1 += x(lfsu1,i)*phi1[i];

        // evaluate gradient of basis functions on reference element
        std::vector<JacobianType> js0(lfsu0.size());
        lfsu0.finiteElement().localBasis().evaluateJacobian(it->position(),js0);
        std::vector<JacobianType> js1(lfsu1.size());
        lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js1);

        // transform gradients from reference element to real element
        const typename EG::Geometry::JacobianInverseTransposed
          jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > gradphi0(lfsu0.size());
        for (size_type i=0; i<lfsu0.size(); i++)
          jac.mv(js0[i][0],gradphi0[i]);
        std::vector<Dune::FieldVector<RF,dim> > gradphi1(lfsu1.size());
        for (size_type i=0; i<lfsu1.size(); i++)
          jac.mv(js1[i][0],gradphi1[i]);

        // compute gradient of u_0, u_1
        Dune::FieldVector<RF,dim> gradu0(0.0);
        for (size_type i=0; i<lfsu0.size(); i++)
          gradu0.axpy(x(lfsu0,i),gradphi0[i]);
        Dune::FieldVector<RF,dim> gradu1(0.0);
        for (size_type i=0; i<lfsu1.size(); i++)
          gradu1.axpy(x(lfsu1,i),gradphi1[i]);

        // integrate both components
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());
        // eq. 0: - d_0 \Delta u_0 - (\lambda*u_0 - u_0^3 - \sigma* u_1 + \kappa) = 0
        for (size_type i=0; i<lfsu0.size(); i++)
          r.accumulate(lfsu0,i,(d_0*(gradu0*gradphi0[i])
				-(lambda*u_0-u_0*u_0*u_0-sigma*u_1+kappa)
				*phi0[i])*factor);
        // eq. 1: - d_1 \Delta u_1 - (u_0 - u_1) = 0
        for (size_type i=0; i<lfsu1.size(); i++)
          r.accumulate(lfsu1,i,(d_1*(gradu1*gradphi1[i])
				-(u_0-u_1)*phi1[i])*factor);
      }
  }

private:
  unsigned int intorder;
  double d_0, d_1, lambda, sigma, kappa;
};
