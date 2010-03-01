// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_COMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_COMPONENTTRANSPORTOP_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/grid/common/genericreferenceelements.hh>

#include<dune/localfunctions/raviartthomas/raviartthomas0q.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/transportccfv.hh>

namespace Dune {
  namespace PDELab {

    //! base class for parameter class
	template<class T, class Imp>
	class ModifiedTransportSpatialParameterInterface 
      : public TransportSpatialParameterInterface<T,Imp>
	{
	public:
	  typedef T Traits;

	  //! saturation at end of big step
	  typename Traits::RangeFieldType 
	  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().snew(e,x);
	  }

	  //! saturation at end of big step
	  typename Traits::RangeFieldType 
	  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().sold(e,x);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};


	/** a local operator for a cell-centered finite folume scheme for
        the transport equation
	  
		\nabla \cdot \{v u - D \nabla u \} = q in \Omega
                                         u = g on \Gamma_D
           \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
                                       outflow on \Gamma_O        

        Modified version for the case

		 d_t (c(x,t)u(x,t)) + \nabla \cdot \{v u - D \nabla u \} = q in \Omega

        where c(x,t) may become zero. We assume that the following holds:

        c(x,t+dt) <= eps  ==>  c(x,t) <= eps

		\tparam TP  parameter class implementing ComponentTransportParameterInterface
	*/
    template<typename TP>
	class ModifiedCCFVSpatialTransportOperator : 
      //      public NumericalJacobianSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      //      public NumericalJacobianBoundary<ModifiedCCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianApplySkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianApplyBoundary<ModifiedCCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianApplyVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
	{
      enum { dim = TP::Traits::GridViewType::dimension };

	public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

	  // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
	  enum { doAlphaVolumePostSkeleton = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume    = true };

      enum { doSkeletonTwoSided = true }; // need to see face from both sides for CFL calculation

      ModifiedCCFVSpatialTransportOperator (TP& tp_) 
		: tp(tp_), zero(1e-7)
	  {
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

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);
        
        // evaluate saturation at end of big step for cell activity
        snew = tp.snew(eg.entity(),inside_local);
        active_cell = snew>zero;
        if (active_cell) 
          active_cell_count++;
        else
          inactive_cell_count++;

        celloutflux = 0.0; // prepare dt computation

//           std::cout << "alpha_volume: time=" << time
//                     << " pos=" << eg.geometry().center()
//                     << " snew=" << snew
//                     << " inactive_cell_count=" << inactive_cell_count
//                     << " active_cell_count=" << active_cell_count
//                     << " active_cell=" << active_cell
//                     << " x=" << x[0]
//                     << std::endl;
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        snew = tp.snew(eg.entity(),inside_local);
        active_cell = snew>zero;
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

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        // activity of neighbor
        bool active_neighbor = tp.snew(*(ig.outside()),outside_local)>zero;

        // convective flux
        RF u_upwind=0.0;
        if (vn>=0) 
          u_upwind = x_s[0]; 
        else 
          {
            u_upwind = x_n[0];
          }
        r_s[0] += (u_upwind*vn)*face_volume;
        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

        // evaluate diffusion coefficients
        typename TP::Traits::RangeFieldType D_inside = tp.sold(*(ig.inside()),inside_local)*tp.D(*(ig.inside()),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.sold(*(ig.outside()),outside_local)*tp.D(*(ig.outside()),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
 
        // diffusive flux
        // note: we do only one-sided evaluation here
        r_s[0] -= (D_avg*(x_n[0]-x_s[0])/distance)*face_volume;
        celloutflux += D_avg*face_volume/distance;
      }


      // jacobian of skeleton term
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_skeleton (const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn, 
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>& 
          outside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        // activity of neighbor
        bool active_neighbor = tp.snew(*(ig.outside()),outside_local)>zero;

        // convective flux
        RF u_upwind=0.0;
        if (vn>=0) 
          {
            mat_ss(0,0) += vn*face_volume;
          }
        else 
          {
            mat_sn(0,0) += vn*face_volume;
          }

        // evaluate diffusion coefficients
        typename TP::Traits::RangeFieldType D_inside = tp.sold(*(ig.inside()),inside_local)*tp.D(*(ig.inside()),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.sold(*(ig.outside()),outside_local)*tp.D(*(ig.outside()),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension> 
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension> 
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();
 
        // diffusive flux
        // note: we do only one-sided evaluation here
        mat_ss(0,0) += D_avg/distance*face_volume;
        mat_sn(0,0) -= D_avg/distance*face_volume;
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

        // no fluxes if cell is not active !
        if (!active_cell) return;
    
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        int bc = tp.bc(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bc==0) // Neumann boundary
          {
            typename TP::Traits::RangeFieldType j = tp.j(*(ig.inside()),face_center_in_element);
            r_s[0] += j*face_volume;
            return;
         }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

        if (bc==2) // Outflow boundary
          {
            r_s[0] += vn*x_s[0]*face_volume;
            return;
          }

        if (bc==1) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            if (vn>=0) g=x_s[0]; else g=tp.g(*(ig.inside()),face_center_in_element);
            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.sold(*(ig.inside()),inside_local)*tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension> 
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension> 
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            r_s[0] += (g*vn - D_inside*(g-x_s[0])/distance)*face_volume;
            celloutflux += D_inside*face_volume/distance;
            return;
          }
      }

      // jacobian of boundary term
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_boundary (const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              LocalMatrix<R>& mat_ss) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;

        // no fluxes if cell is not active !
        if (!active_cell) return;
    
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>& 
          face_local = Dune::GenericReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        int bc = tp.bc(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bc==0) // Neumann boundary
          {
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        if (bc==2) // Outflow boundary
          {
            mat_ss(0,0) += vn*face_volume;
            return;
          }

        if (bc==1) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            if (vn>=0) 
              mat_ss(0,0) += vn*face_volume;
            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::GenericReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.sold(*(ig.inside()),inside_local)*tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension> 
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension> 
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            mat_ss(0,0) += D_inside/distance*face_volume;
            return;
          }
      }

      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        if (!first_stage) return; // time step calculation is only done in first stage

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = tp.c(eg.entity(),inside_local)*eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(celloutflux+1E-40);

//         if (active_cell)
//           std::cout << "A: time=" << time
//                     << " pos=" << eg.geometry().center()
//                     << " snew=" << snew
//                     << " capacityvol=" << cellcapacity
//                     << " outflux=" << celloutflux
//                     << " residual=" << r[0]
//                     << " celldt=" << celldt
//                     << " x=" << x[0]
//                     << std::endl;

        if (active_cell)
          dtmin = std::min(dtmin,celldt);
      }


 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
		// domain and range field type
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSV::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        // no sources if cell is not active !
        if (!active_cell) return;
    
        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate source term
        typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local);

        r[0] -= q*eg.geometry().volume();
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt, 
                    int stages)
      {
        tp.preStep(time,dt,stages);
      }
      
      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
//         std::cout << "preStage on transport operator called" << std::endl;
        if (r==1)
          {
            first_stage = true;
            dtmin = 1E100;
            active_cell_count = 0;
            inactive_cell_count = 0;
          }
        else first_stage = false;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
//         std::cout << "active cells: " << active_cell_count << " dtmin: " << dtmin << std::endl;
//        std::cout << "inactive cells: " << inactive_cell_count << " dtmin: " << dtmin << std::endl;
        return dtmin;
      }

	private:

      typename TP::Traits::RangeFieldType regularization (typename TP::Traits::RangeFieldType x)
      {
        const typename TP::Traits::RangeFieldType regeps = 1e-20;
        const typename TP::Traits::RangeFieldType regeps2 = regeps*regeps;
        if (x<-regeps) return -1;
        if (x<0.0) return (x+regeps)*(x+regeps)/regeps2-1.0;
        if (x<regeps) return 1.0-(regeps-x)*(regeps-x)/regeps2;
        return +1.0;
      }

	  TP& tp;
      bool first_stage;
      mutable bool active_cell;
      mutable typename TP::Traits::RangeFieldType snew;
      mutable int active_cell_count;
      mutable int inactive_cell_count;
      typename TP::Traits::RangeFieldType time;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      mutable typename TP::Traits::RangeFieldType celloutflux;
      typename TP::Traits::RangeFieldType zero;
	};


    //! base class for parameter class
	template<class T, class Imp>
	class ModifiedTransportTemporalParameterInterface 
      : public TransportTemporalParameterInterface<T,Imp>
	{
	public:
	  typedef T Traits;

	  //! saturation at end of big step
	  typename Traits::RangeFieldType 
	  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
	  {
		return asImp().snew(e,x);
	  }

	private:
	  Imp& asImp () {return static_cast<Imp &> (*this);}
	  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
	};

   /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x,t) uv dx
     * \f}
     *
     * version where c(x,t) may become zero.
     */
    template<class TP>
	class ModifiedCCFVTemporalOperator 
      : public NumericalJacobianApplyVolume<ModifiedCCFVTemporalOperator<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      ModifiedCCFVTemporalOperator (TP& tp_) 
		: tp(tp_), zero(1e-7)
	  {
	  }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
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

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local);
        bool active_cell = snew>zero;

        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);

        if (active_cell)
          {
            // residual contribution
            r[0] += c*x[0]*eg.geometry().volume();
          }

//         if (active_cell)
//           std::cout << "M: time=" << time
//                   << " pos=" << eg.geometry().center()
//                   << " snew=" << snew
//                   << " capacityvol=" << c*eg.geometry().volume()
//                   << " residual=" << r[0]
//                   << std::endl;
	  }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, 
                            LocalMatrix<R>& mat) const
      {
		// domain and range field type
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

        // cell center
        const Dune::FieldVector<DF,dim>& 
          inside_local = Dune::GenericReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local);
        bool active_cell = snew>zero;

        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);
        
        // residual contribution
        if (active_cell) 
          mat(0,0) += c*eg.geometry().volume();
        else
          mat(0,0) += 1.0; // rhs should be zero so we get zero update

//         if (active_cell)
//           std::cout << "J: time=" << time
//                   << " pos=" << eg.geometry().center()
//                   << " snew=" << snew
//                   << " mat=" << mat(0,0)
//                   << std::endl;
      }

	private:
	  TP& tp;
      typename TP::Traits::RangeFieldType time;
      typename TP::Traits::RangeFieldType zero;
	};

  }
}

#endif
