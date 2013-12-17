#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with cell-centered finite volumes on axiparallel, structured grids
 *
 * \tparam B a function indicating the type of boundary condition
 * \tparam G a function for the values of the Dirichlet boundary condition
 */
class Example04LocalOperator :    // implement jacobian evaluation in base classes
  public Dune::PDELab::NumericalJacobianApplyVolume<Example04LocalOperator>,
  public Dune::PDELab::NumericalJacobianVolume<Example04LocalOperator>,
  public Dune::PDELab::NumericalJacobianApplySkeleton<Example04LocalOperator>,
  public Dune::PDELab::NumericalJacobianSkeleton<Example04LocalOperator>,
  public Dune::PDELab::NumericalJacobianApplyBoundary<Example04LocalOperator>,
  public Dune::PDELab::NumericalJacobianBoundary<Example04LocalOperator>,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton 
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume  = true };
  enum { doAlphaSkeleton  = true };                             // assemble skeleton term
  enum { doAlphaBoundary  = true };

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
		     R& r) const
  {
    // range field type
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    // evaluate reaction term
    RF a = 0.0;
    RF f = 0.0;

    r.accumulate(lfsu,0,(a*x(lfsu,0)-f)*eg.geometry().volume());
  }

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE! 
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig, 
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                       R& r_s, R& r_n) const
  {
    // domain and range field type
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    const int dim = IG::dimension;

    // distance between cell centers in global coordinates
    Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().center();
    Dune::FieldVector<DF,dim> outside_global = ig.outside()->geometry().center();
    inside_global -= outside_global;
    RF distance = inside_global.two_norm();

    // face geometry
    RF face_volume = ig.geometry().volume();
 
    // diffusive flux for both sides
    r_s.accumulate(lfsu_s,0,-(x_n(lfsu_n,0)-x_s(lfsu_s,0))*face_volume/distance);
    r_n.accumulate(lfsu_n,0,(x_n(lfsu_n,0)-x_s(lfsu_s,0))*face_volume/distance);
  }

  // skeleton integral depending on test and ansatz functions
  // Here Dirichlet and Neumann boundary conditions are evaluated
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, 
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       R& r_s) const
  {
    // domain and range field type
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    const int dim = IG::dimension;
    
    // face geometry
    Dune::FieldVector<DF,dim> face_center = ig.geometry().center();
    RF face_volume = ig.geometry().volume();

    // evaluate boundary condition type
    int b;
    if (face_center[0]>1.0-1e-6)
      b = 0; // Neumann
    else
      b = 1; // Dirichlet

    // do things depending on boundary condition type
    if (b==0) // Neumann boundary
      {
        RF j; if (face_center[1]<0.5) j = 1.0; else j = -1.0;
        r_s.accumulate(lfsu_s,0,j*face_volume);
        return;
      }

    if (b==1) // Dirichlet boundary
      {
        RF g;
        if (face_center[0]<1E-6 && face_center[1]>0.25 && face_center[1]<0.5)
          g = 1.0;
        else
          g = 0.0;
        Dune::FieldVector<DF,dim> inside_global = ig.inside()->geometry().center();
        inside_global -= face_center;
        RF distance = inside_global.two_norm();
        r_s.accumulate(lfsu_s,0,-(g-x_s(lfsu_s,0))*face_volume/distance);
        return;
      }
  }
};
