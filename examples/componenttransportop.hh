// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_COMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_COMPONENTTRANSPORTOP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>

#include<dune/finiteelements/rt0q.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {

	//! traits class for two phase parameter class
	template<typename GV, typename RF>
	struct ComponentTransportParameterTraits
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

    //! base class for parameter class
	template<class T, class Imp>
	class ComponentTransportParameterInterface
	{
	public:
	  typedef T Traits;

	  //! porosity
	  typename Traits::RangeFieldType 
	  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().phi(e,x);
	  }

	  //! saturation at new time level
	  typename Traits::RangeFieldType 
	  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().snew(e,x);
	  }

	  //! saturation at old time level
	  typename Traits::RangeFieldType 
	  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().sold(e,x);
	  }

	  //! velocity field (at new time level)
	  typename Traits::RangeType
	  u (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().u(e,x);
	  }

	  //! Dirichlet boundary condition on inflow
	  typename Traits::RangeFieldType 
	  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, 
         typename Traits::RangeFieldType time) const
	  {
		return asImp().g(is,x,time);
	  }

	  //! source term
	  typename Traits::RangeFieldType 
	  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x, 
         typename Traits::RangeFieldType time) const
	  {
		return asImp().q(e,x,time);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

	/** a local operator for an explicit finite volume scheme
		to compute transport of a dissolved component in a phase
	  
		\partial_t (\phi s c) + \nabla \cdot \{c u} - q = 0 in \Omega \times I 

		- c is the unknown concentration
		- Implements explicit Euler scheme with full upwinding. 
		- Saturation s may become zero
		- Provides CFL condition, actual update is done externally.
		- We use only residual evaluation, no matrix can be assembled

		TP : parameter class implementing ComponentTransportParameterInterface
		V  : Vector holding last time step
	*/
    template<typename TP, typename V>
	class EEComponentTransportOperator : public LocalOperatorDefaultFlags
	{
      enum { dim = TP::Traits::GridViewType::dimension };

	public:
	  // residual assembly flags
      enum { doAlphaVolume    = true };
	  enum { doAlphaVolumePostSkeleton = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doSkeletonTwoSided = true }; // need to see face from both sides for CFL calculation

      EEComponentTransportOperator (const TP& tp_, V& scaling_) 
		: tp(tp_), scaling(scaling_)
	  {
		time = 0;
        zero_saturation = 1E-8;
        min_snew, min_sold, min_local_flux_sum = 0;
        cellno = -1;
	  }

	  // set time where operator is to be evaluated (i.e. end of the time intervall)
	  void set_time (typename TP::Traits::RangeFieldType time_)
	  {
		time = time_;
	  } 

	  // get maximum allowable time step computed while assembling residual
	  typename TP::Traits::RangeFieldType get_max_timestep () const
	  {
            std::cout << "minimum time step info: cell No. " << cellno
                      << " sold=" << min_sold
                      << " snew=" << min_snew
                      << " local_flux_sum=" << min_local_flux_sum << std::endl;
		return max_timestep;
	  } 

	  // get number of cells where saturation is zero
	  int get_empty_cells () const
	  {
		return emptycells;
	  } 

	  // to be called before assembling residual
	  void initialize_timestep ()
	  {
		max_timestep = 1E100;
        emptycells = 0;
	  } 

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // initialize local time step calculation
        local_flux_sum = 0;

        // cell geometry
        const Dune::FieldVector<DF,dim>& 
          cell_center_local = Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        cell_volume = eg.geometry().integrationElement(cell_center_local)
          *Dune::GenericReferenceElements<DF,dim>::general(eg.geometry().type()).volume();

        // porosity and saturation
		porosity = tp.phi(eg.entity(),cell_center_local);
		snew = tp.snew(eg.entity(),cell_center_local);
		sold = tp.sold(eg.entity(),cell_center_local);

        if (snew>zero_saturation && sold>zero_saturation)
          activecell = true;
        else
          activecell = false;

//         std::cout << "cell No. " << lfsu.globalIndex(0) 
//                   << " porosity=" << porosity 
//                   << " cell_volume=" << cell_volume 
//                   << std::endl;

		// scale value of concentration at old time level
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;
        RF q = tp.q(eg.entity(),cell_center_local,time); // can not eval at new time because timestep is not known
        if (activecell)
          {
            B::access(scaling,lfsu.globalIndex(0)) = sold/snew;
            r[0] -= q/(porosity*snew);
          }
        else
          {
            //           std::cout << "Bingo ! cell " << lfsu.globalIndex(0) << std::endl;
            B::access(scaling,lfsu.globalIndex(0)) = 0;
            emptycells++;
          }
      }

      // post skeleton: compute time step allowable for cell
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        RF local_timestep = sold/(std::abs(local_flux_sum)+1E-30);
        if (activecell && local_timestep<max_timestep) 
          {
            max_timestep = local_timestep;
            min_sold = sold;
            min_snew = snew;
            min_local_flux_sum = local_flux_sum;
            cellno = lfsu.globalIndex(0);
          }

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
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();

        // face center in element coordinates
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType u(tp.u(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF wn = u*ig.unitOuterNormal(face_local);

//         std::cout << "cell No. " << lfsu_s.globalIndex(0) 
//                   << " u=" << u 
//                   << " n=" << ig.unitOuterNormal(face_local)
//                   << " wn=" << wn 
//                   << std::endl;

        // upwind
        RF c_upwind;
        if (wn>=0)
          c_upwind = x_s[0];
        else
          c_upwind = x_n[0];

        // residual updates (use asymetric flux evaluation, may result in nonconservation of stuff
        // note. we do only one-sided evaluation here
        if (activecell)
          {
            r_s[0] += c_upwind*wn*face_volume/(cell_volume*porosity*snew);
            if (wn>0) local_flux_sum += wn*face_volume/(cell_volume*porosity);
          }
      }

	  // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary (const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
	  {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;

        if (!activecell) return; // empty cell; nothing to do
    
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).volume();

        // face center in element coordinates
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType u(tp.u(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF wn = u*ig.unitOuterNormal(face_local);

        // evaluate boundary condition
        RF c_boundary = tp.g(ig.intersection(),face_local,time); // this is incorrect, but we do not know the timestep!!

        // residual updates (use asymetric flux evaluation, may result in nonconservation of stuff
        // note. we do only one-sided evaluation here
        if (wn>=0) 
          { 
            // outflow boundary
            r_s[0] += x_s[0]*wn*face_volume/(cell_volume*porosity*snew);
            local_flux_sum += wn*face_volume/(cell_volume*porosity);
          }
        else
          {
            // inflow boundary
            r_s[0] += c_boundary*wn*face_volume/(cell_volume*porosity*snew);
          }
      }

	private:
	  const TP& tp;
	  V& scaling;
	  typename TP::Traits::RangeFieldType time;
	  mutable typename TP::Traits::RangeFieldType max_timestep; 
	  mutable typename TP::Traits::RangeFieldType local_flux_sum;
      mutable typename TP::Traits::RangeFieldType cell_volume;
      mutable typename TP::Traits::RangeFieldType porosity;
      mutable typename TP::Traits::RangeFieldType snew,sold;
      typename TP::Traits::RangeFieldType zero_saturation;
      mutable int emptycells;
      mutable bool activecell;
      mutable int cellno;
      mutable typename TP::Traits::RangeFieldType min_snew, min_sold, min_local_flux_sum;
	};

  }
}

#endif
