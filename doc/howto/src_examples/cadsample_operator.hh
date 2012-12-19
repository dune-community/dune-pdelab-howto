#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

/** a local operator for solving the equation
 *
 *   - \Delta m(x) u + a*u   = f   in \Omega
 *                       u   = g   on \Gamma_D\subseteq\partial\Omega
 *   - \nabla m(x) u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 * Note the source term is 0 here.
 *
 * \tparam M a function indicating the diffusion coefficient
 * \tparam B a function indicating the type of boundary condition
 * \tparam J a function indicating the source term
 */
template<typename M, typename B, typename J>                            // NEW
class CADLocalOperator :

  public Dune::PDELab::NumericalJacobianApplyVolume<CADLocalOperator<M, B, J> >,
  public Dune::PDELab::NumericalJacobianVolume<CADLocalOperator<M, B, J> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<CADLocalOperator<M, B, J> >,
  public Dune::PDELab::NumericalJacobianBoundary<CADLocalOperator<M, B, J> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{

public:

  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };

  // constructor parametrized by material and boundary classes          // NEW
  CADLocalOperator (const M& m_, const B& b_, const J& j_, unsigned int intorder_=2)
    : m(m_), b(b_), j(j_), intorder(intorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // extract some types
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator
           it=rule.begin(); it!=rule.end(); ++it)
    {
      // evaluate basis functions on reference element
      std::vector<RangeType> phi(lfsu.size());
      lfsu.localFiniteElement().localBasis().evaluateFunction(it->position(),phi);

      // compute u at integration point
      RF u=0.0;
      for (size_type i=0; i<lfsu.size(); i++)
        u += x[i]*phi[i];

      // evaluate gradient of basis functions on reference element
      std::vector<JacobianType> js(lfsu.size());
      lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

      // transform gradients from reference element to real element
      const typename EG::Geometry::JacobianInverseTransposed
        jac = eg.geometry().jacobianInverseTransposed(it->position());
      std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
      for (size_type i=0; i<lfsu.size(); i++)
        jac.mv(js[i][0],gradphi[i]);

      // compute gradient of u
      Dune::FieldVector<RF,dim> gradu(0.0);
      for (size_type i=0; i<lfsu.size(); i++)
        gradu.axpy(x[i],gradphi[i]);

      // evaluate parameters;
      Dune::FieldVector<RF,dim>
        globalpos = eg.geometry().global(it->position());
      RF f = 0;
      RF a = 0;
      typename M::Traits::RangeType y;
      m.evaluate(eg.entity(), it->position(), y);                       // NEW
      gradu *= (double) y;

      // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
      RF factor = it->weight()*eg.geometry().integrationElement(it->position());
      for (size_type i=0; i<lfsu.size(); i++)
        r[i] += ( gradu*gradphi[i] + a*u*phi[i] - f*phi[i] )*factor;
    }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // some types
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSV::Traits::LocalFiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSV::Traits::SizeType size_type;

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate boundary condition type
        typename B::Traits::RangeType bctype;
        b.evaluate(ig,it->position(),bctype);

        // skip rest if we are on Dirichlet boundary
        if (bctype>0) continue;

        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate basis functions at integration point
        std::vector<RangeType> phi(lfsv_s.size());
        lfsu_s.localFiniteElement().localBasis().evaluateFunction(local,phi);

        // evaluate u (e.g. flux may depend on u)
        RF u=0.0;
        for (size_type i=0; i<lfsu_s.size(); i++)
          u += x_s[i]*phi[i];

        // evaluate flux boundary condition
        typename J::Traits::RangeType y;                                // NEW
        j.evaluate(ig, it->position(), y);

        // integrate j
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s[i] += y * phi[i] * factor;
      }
  }

private:

  const M& m;
  const B& b;
  const J& j;
  unsigned int intorder;
};
