// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TWOPHASEOP_HH
#define DUNE_PDELAB_TWOPHASEOP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {

	//! traits class for two phase parameter class
	template<typename GV, typename RF>
	struct TwoPhaseParameterTraits
	{
	  //! \brief the grid view
	  typedef GV GridViewType;

	  //! \brief Enum for domain dimension
	  enum { 
		//! \brief dimension of the domain
		dimDomain = GV::dimension
	  }; 

	  //! \brief Export type for domain field
	  typedef typename GV::Grid::ctype DomainFieldType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

	  //! \brief domain type
	  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

	  //! \brief Export type for range field
	  typedef RF RangeFieldType;

	  //! \brief range type
	  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

	  //! grid types
	  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
	  typedef typename GV::Intersection IntersectionType;
	};

	template<class T, class Imp>
	class TwoPhaseParameterInterface
	{
	public:
	  typedef T Traits;

	  //! porosity
	  typename Traits::RangeFieldType 
	  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().phi(e,x);
	  }

	  //! capillary pressure function
	  typename Traits::RangeFieldType 
	  pc (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
		   typename Traits::RangeFieldType s_l) const
	  {
		return asImp().pc(e,x,s_l);
	  }
	  
	  //! inverse capillary pressure function
	  typename Traits::RangeFieldType 
	  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
		   typename Traits::RangeFieldType pc) const
	  {
		return asImp().s_l(e,x,pc);
	  }
	  
	  //! liquid phase relative permeability
	  typename Traits::RangeFieldType 
	  kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType s_l) const
	  {
		return asImp().kr_l(e,x,s_l);
	  }

	  //! gas phase relative permeability
	  typename Traits::RangeFieldType 
	  kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType s_g) const
	  {
		return asImp().kr_g(e,x,s_g);
	  }

	  //! liquid phase viscosity
	  typename Traits::RangeFieldType 
	  mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType p_l) const
	  {
		return asImp().mu_l(e,x,p_l);
	  }

	  //! gas phase viscosity
	  typename Traits::RangeFieldType 
	  mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType p_g) const
	  {
		return asImp().mu_l(e,x,p_g);
	  }
	  
	  //! absolute permeability (scalar!)
	  typename Traits::RangeFieldType 
	  k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().k_abs(e,x);
	  }

	  //! gravity vector
	  const typename Traits::RangeType& gravity () const
	  {
		return asImp().gravity();
	  }

	  //! liquid phase molar density
	  template<typename E>
	  typename Traits::RangeFieldType 
	  nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType p_l) const
	  {
		return asImp().nu_l(e,x,p_l);
	  }

	  //! gas phase molar density
	  typename Traits::RangeFieldType 
	  nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			typename Traits::RangeFieldType p_g) const
	  {
		return asImp().nu_g(e,x,p_g);
	  }

	  //! liquid phase mass density
	  typename Traits::RangeFieldType 
	  rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			 typename Traits::RangeFieldType p_l) const
	  {
		return asImp().rho_l(e,x,p_l);
	  }

	  //! gas phase mass density
	  typename Traits::RangeFieldType 
	  rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
			 typename Traits::RangeFieldType p_g) const
	  {
		return asImp().rho_g(e,x,p_g);
	  }
	  
	  //! liquid phase boundary condition type
	  int
	  bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().bc_l(is,x);
	  }

	  //! gas phase boundary condition type
	  int
	  bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().bc_g(is,x);
	  }

	  //! liquid phase Dirichlet boundary condition
	  typename Traits::RangeFieldType 
	  g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().g_l(is,x);
	  }

	  //! gas phase Dirichlet boundary condition
	  typename Traits::RangeFieldType 
	  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().g_g(is,x);
	  }

	  //! liquid phase Neumann boundary condition
	  typename Traits::RangeFieldType 
	  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().j_l(is,x);
	  }

	  //! gas phase Neumann boundary condition
	  typename Traits::RangeFieldType 
	  j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
	  {
		return asImp().j_g(is,x);
	  }

	  //! liquid phase source term
	  typename Traits::RangeFieldType 
	  q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
		   typename Traits::RangeFieldType time) const
	  {
		return asImp().q_l(e,x,time);
	  }

	  //! gas phase source term
	  typename Traits::RangeFieldType 
	  q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
		   typename Traits::RangeFieldType time) const
	  {
		return asImp().q_g(e,x,time);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


	// a local operator for solving the two-phase flow in pressure-pressure formulation
	// with two-point flux approximation
    // TP : parameter class, see above
	// V  : Vector holding last time step
    template<typename TP, typename V>
	class TwoPhaseTwoPointFluxOperator : public NumericalJacobianVolume<TwoPhaseTwoPointFluxOperator<TP,V> >,
										 public NumericalJacobianApplyVolume<TwoPhaseTwoPointFluxOperator<TP,V> >,

										 public NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP,V> >,
										 public NumericalJacobianApplySkeleton<TwoPhaseTwoPointFluxOperator<TP,V> >,

										 public NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP,V> >,
										 public NumericalJacobianApplyBoundary<TwoPhaseTwoPointFluxOperator<TP,V> >,

										 public FullSkeletonPattern, 
										 public FullVolumePattern,
                                         public LocalOperatorDefaultFlags
	{
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = 0 };
      enum { gas = 1 };

	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

	  // residual assembly flags
      enum { doAlphaVolume    = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      enum { doLambdaBoundary = true };

      TwoPhaseTwoPointFluxOperator (const TP& tp_, const V& pold_) : tp(tp_), pold(pold_) {}

	  // set time where operator is to be evaluated (i.e. end of the time intervall)
	  void set_time (typename TP::Traits::RangeFieldType time_)
	  {
		time = time_;
	  } 

	  // set size of the time step in implicit Euler
	  void set_timestep (typename TP::Traits::RangeFieldType timestep_)
	  {
		timestep = timestep_;
	  } 

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;
 
		// domain and range field type
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometry
        const Dune::FieldVector<DF,dim>& 
          cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().integrationElement(cell_center_local)
          *Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).volume();

		RF phi = tp.phi(eg.entity(),cell_center_local);
		RF s_l = tp.s_l(eg.entity(),cell_center_local,x[gas]-x[liquid]);

		r[liquid] += phi * s_l * tp.nu_l(eg.entity(),cell_center_local,x[liquid]) * cell_volume;
		r[gas]    += phi * (1-s_l) * tp.nu_g(eg.entity(),cell_center_local,x[gas]) * cell_volume;
	  }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<liquid>::Type PLSpace;
		const PLSpace& plspace = lfsv.template getChild<liquid>();
        typedef typename LFSV::template Child<gas>::Type PGSpace;
		const PGSpace& pgspace = lfsv.template getChild<gas>();

		// domain and range field type
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometry
        const Dune::FieldVector<DF,dim>& 
          cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().integrationElement(cell_center_local)
          *Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).volume();

		// contribution from last time step
		typedef typename LFSV::Traits::GridFunctionSpaceType::Traits::BackendType B;
		RF p_l_old = B::const_access(pold,lfsv.globalIndex(liquid));
		RF p_g_old = B::const_access(pold,lfsv.globalIndex(gas));
		RF phi = tp.phi(eg.entity(),cell_center_local);
		RF s_l = tp.s_l(eg.entity(),cell_center_local,p_g_old-p_l_old);
		r[liquid] -= phi * s_l * tp.nu_l(eg.entity(),cell_center_local,p_l_old) * cell_volume;
		r[gas]    -= phi * (1-s_l) * tp.nu_g(eg.entity(),cell_center_local,p_g_old) * cell_volume;

		// contribution from source term
		r[liquid] -= tp.q_l(eg.entity(),cell_center_local,time) * timestep * cell_volume;
		r[gas]    -= tp.q_g(eg.entity(),cell_center_local,time) * timestep * cell_volume;
	  }

	  // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n) const
	  {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;
 
		// domain and range field type
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // cell geometries
        const Dune::FieldVector<DF,dim>& 
          inside_cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>& 
          outside_cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        Dune::FieldVector<DF,IG::dimension> 
          inside_cell_center_global = ig.inside()->geometry().global(inside_cell_center_local);
        Dune::FieldVector<DF,IG::dimension> 
          outside_cell_center_global = ig.outside()->geometry().global(outside_cell_center_local);

        // distance of cell centers
        Dune::FieldVector<DF,dim> d(outside_cell_center_global);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);
        RF k_abs_outside = tp.k_abs(*(ig.outside()),outside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]);
        RF rho_l_outside = tp.rho_l(*(ig.outside()),outside_cell_center_local,x_n[liquid]);
        RF w_l = (x_s[liquid]-x_n[liquid])/distance + aavg(rho_l_inside,rho_l_outside)*gn; // determines direction
        RF pc_upwind, s_l_upwind, s_g_upwind;
        if (w_l>=0) // upwind capillary pressure on face
          {
            pc_upwind = x_s[gas]-x_s[liquid];
            s_l_upwind = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
          }
        else
          {
            pc_upwind = x_n[gas]-x_n[liquid];
            s_l_upwind = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
          }
        s_g_upwind = 1-s_l_upwind;
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]);
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n[liquid]);
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);
        RF nu_l = aavg(tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]), 
                       tp.nu_l(*(ig.outside()),outside_cell_center_local,x_n[liquid]));
//         std::cout << "   timestep = " << timestep << std::endl;
//         std::cout << "       nu_l = " << nu_l << std::endl;
//         std::cout << "    sigma_l = " << sigma_l << std::endl;
//         std::cout << "        w_l = " << w_l << std::endl;
//         std::cout << "face_volume = " << face_volume << std::endl;
//         std::cout << std::endl;

        r_s[liquid] += timestep * nu_l * sigma_l * w_l * face_volume;
        r_n[liquid] -= timestep * nu_l * sigma_l * w_l * face_volume;

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s[gas]);
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n[gas]);
        RF w_g = (x_s[gas]-x_n[gas])/distance + aavg(rho_g_inside,rho_g_outside)*gn; // determines direction
        if (w_l*w_g<0) // new evaluation necessary only if signs differ
          {
            if (w_g>=0) // upwind capillary pressure on face
              {
                pc_upwind = x_s[gas]-x_s[liquid];
                s_l_upwind = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
              }
            else
              {
                pc_upwind = x_n[gas]-x_n[liquid];
                s_l_upwind = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
              }
            s_g_upwind = 1-s_l_upwind;
         }
        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_upwind)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s[gas]);
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_upwind)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n[gas]);
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);
        RF nu_g = aavg(tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s[gas]), 
                       tp.nu_g(*(ig.outside()),outside_cell_center_local,x_n[gas]));

        r_s[gas] += timestep * nu_g * sigma_g * w_g * face_volume;
        r_n[gas] -= timestep * nu_g * sigma_g * w_g * face_volume;
	  }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
	  {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;
 
		// domain and range field type
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // face geometry
        const Dune::FieldVector<DF,dim-1>& 
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::GenericReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

        // cell geometry
        const Dune::FieldVector<DF,dim>& 
          inside_cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        Dune::FieldVector<DF,dim> 
          inside_cell_center_global = ig.inside()->geometry().global(inside_cell_center_local);

        // distance of cell center to boundary
        Dune::FieldVector<DF,dim> d = ig.geometry().global(face_local);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase Dirichlet boundary
        if (bc_l==1) 
          {
            RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]);
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s[liquid]-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s[gas]-x_s[liquid]);
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]);
            RF sigma_l = lambda_l_inside*k_abs_inside;
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s[liquid]);
            r_s[liquid] += timestep * nu_l * sigma_l * w_l * face_volume;
          }

        // gas phase Dirichlet boundary
        if (bc_g==1) 
          {
            RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s[gas]);
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s[gas]-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s[gas]-x_s[liquid]);
            RF s_g = 1-s_l;
            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s[gas]);
            RF sigma_g = lambda_g_inside*k_abs_inside;
            RF nu_g = tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s[gas]);
            r_s[gas] += timestep * nu_g * sigma_g * w_g * face_volume;
          }
      }

      // boundary integral independent of ansatz functions
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;
 
		// domain and range field type
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PLSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // face geometry
        const Dune::FieldVector<DF,dim-1>& 
          face_local = Dune::GenericReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::GenericReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=0 && bc_g!=0) return; // no Neumann boundary conditions

        // liquid phase Neumann boundary
        if (bc_l==0) 
          {
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s[liquid] += timestep * j_l * face_volume;
          }
 
        // gas phase Neumann boundary
        if (bc_g==0) 
          {
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s[gas] += timestep * j_g * face_volume;
          }
      }

    private:
      const TP& tp;  // two phase parameter class
	  const V& pold; // last time step
	  typename TP::Traits::RangeFieldType time;
	  typename TP::Traits::RangeFieldType timestep;

      template<typename T>
      T aavg (T a, T b) const
      {
        return 0.5*(a+b);
      }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-30;
        return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
      }

	};
  }
}

#endif
