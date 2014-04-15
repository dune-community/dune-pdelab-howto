#ifndef TWOPHASEMFD_HH
#define TWOPHASEMFD_HH

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/mfdcommon.hh>

/* TODO
 * consider molar densities
 */

namespace Dune
{
    namespace PDELab
    {
        template<class Data, class V,
                 class WBuilder = MimeticBrezziW<typename Data::Traits::DomainFieldType,
                                                 Data::Traits::dimDomain> >
        class TwoPhaseMFD
            : public NumericalJacobianVolume<TwoPhaseMFD<Data,V,WBuilder> >
            , public NumericalJacobianApplyVolume<TwoPhaseMFD<Data,V,WBuilder> >

            , public NumericalJacobianSkeleton<TwoPhaseMFD<Data,V,WBuilder> >
            , public NumericalJacobianApplySkeleton<TwoPhaseMFD<Data,V,WBuilder> >

            , public NumericalJacobianBoundary<TwoPhaseMFD<Data,V,WBuilder> >
            , public NumericalJacobianApplyBoundary<TwoPhaseMFD<Data,V,WBuilder> >

            , public FullVolumePattern
            , public FullSkeletonPattern
            , public LocalOperatorDefaultFlags
        {
            static const unsigned int dim = Data::Traits::dimDomain;
            typedef typename Data::Traits::DomainFieldType ctype;
            typedef typename Data::Traits::RangeFieldType rtype;

        public:
            // pattern assembly flags
            enum { doPatternVolume = true };
            enum { doPatternSkeleton = true };

            // residual assembly flags
            enum { doAlphaVolume = true };
            enum { doSkeletonTwoSided = true };
            enum { doAlphaSkeleton = true };
            enum { doAlphaBoundary = true };

            enum { doLambdaVolume = true };
            enum { doLambdaBoundary = true };

            // numbering of the phases
            enum { nonwetting, wetting };

            TwoPhaseMFD(const Data& data_, const V& x_old_, const WBuilder& wbuilder_ = WBuilder())
                : data(data_), x_old(x_old_), wbuilder(wbuilder_), init_mode(false), cell_cache_id(-1)
            {}

            // set time where operator is to be evaluated (i.e. end of the time interval)
            void set_time(rtype time_)
            {
		time = time_;
            }

            // set size of the time step in implicit Euler
            void set_timestep(rtype timestep_)
            {
		timestep = timestep_;
            }

            void set_init_mode(bool init_mode_)
            {
                init_mode = init_mode_;
            }

            // pattern *********************************************

            // define sparsity pattern connecting self and neighbor dofs
            // template<typename LFSU, typename LFSV>
            // void pattern_skeleton (const LFSU& lfsu_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const LFSV& lfsv_n,
            //                        LocalSparsityPattern& pattern_sn, LocalSparsityPattern& pattern_ns) const
            // {
            //     // extract subspaces
            //     typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
            //     typedef typename NonWettingPhase::template Child<0>::Type NonWettingCellUnknowns;
            //     typedef typename LFSV::template Child<wetting>::Type WettingPhase;
            //     typedef typename WettingPhase::template Child<0>::Type WettingCellUnknowns;
            //     const CellUnknowns& cell_space_s = lfsu_s.template getChild<0>();
            //     const CellUnknowns& cell_space_n = lfsu_n.template getChild<0>();
            //     typedef typename LFSV::template Child<1>::Type FaceUnknowns;
            //     const FaceUnknowns& face_space_s = lfsu_s.template getChild<1>();

            //     // add links between all dofs of current cell to cell center dof of neighbor cell
            //     // these links are required for upwinding only
            //     pattern_sn.push_back(SparsityLink(cell_space_s.localIndex(0), cell_space_n.localIndex(0)));
            //     for (unsigned int i = 0; i < face_sapce_s.size(); ++i)
            //     {
            //         pattern_sn.push_back(SparsityLink(face_space_s.localIndex(i), cell_space_n.localIndex(0)));
            //         pattern_ns.push_back(SparsityLink(cell_space_n.localIndex(0), face_space_s.localIndex(i)));
            //     }
            // }

            // alpha ***********************************************

            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
                const NonWettingPhase& n_space = lfsv.template getChild<nonwetting>();
                typedef typename NonWettingPhase::template Child<0>::Type NonWettingCellUnknowns;
                const NonWettingCellUnknowns& n_cell_space = n_space.template getChild<0>();
                typedef typename LFSV::template Child<wetting>::Type WettingPhase;
                const WettingPhase& w_space = lfsv.template getChild<wetting>();
                typedef typename WettingPhase::template Child<0>::Type WettingCellUnknowns;
                const WettingCellUnknowns& w_cell_space = w_space.template getChild<0>();

                // get cell center in local coordinates
                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = GenericReferenceElements<ctype,dim>::general(gt).position(0,0);

                if (cell_cache_id != n_cell_space.globalIndex(0))
                {
                    cell_cache_id = n_cell_space.globalIndex(0);

                    cell.init(eg.entity());

                    typedef typename EG::Entity::LeafIntersectionIterator IntersectionIterator;
                    IntersectionIterator isend = eg.entity().ileafend();
                    for (IntersectionIterator is = eg.entity().ileafbegin(); is != isend; ++is)
                        cell.add_face(*is);

                    // get absolute permeability for current cell
                    const typename Data::Traits::PermTensorType K = data.k_abs(eg.entity(), localcenter);

                    // cache effective gravity vector for current cell
                    K.mv(data.gravity(), K_gravity);

                    // build matrix W
                    wbuilder.build_W(cell, K, W);

                    // cache porosity
                    phi = data.phi(eg.entity(), localcenter);
                }

                // extract wetting phase saturation
                rtype S_w = x(w_cell_space,0);

                if (!init_mode)
                {
                    // storage term (implicit Euler time discretization)
		  r.accumulate(n_cell_space,0,cell.volume * phi * (1.0 - S_w));
		  r.accumulate(w_cell_space,0,cell.volume * phi * S_w);
                }
                else
                {
                    // trivial equations for the cell centres when in initialisation mode
		  r.accumulate(n_cell_space,0,x(n_cell_space,0));
		  r.accumulate(w_cell_space,0,S_w);
                }
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_skeleton(const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                R& r_s, R& r_n) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
                const NonWettingPhase& n_space_s = lfsv_s.template getChild<nonwetting>();
                const NonWettingPhase& n_space_n = lfsv_n.template getChild<nonwetting>();
                typedef typename NonWettingPhase::template Child<0>::Type NonWettingCellUnknowns;
                const NonWettingCellUnknowns& n_cell_space_s = n_space_s.template getChild<0>();
                const NonWettingCellUnknowns& n_cell_space_n = n_space_n.template getChild<0>();
                typedef typename NonWettingPhase::template Child<1>::Type NonWettingFaceUnknowns;
                const NonWettingFaceUnknowns& n_face_space_s = n_space_s.template getChild<1>();

                typedef typename LFSV::template Child<wetting>::Type WettingPhase;
                const WettingPhase& w_space_s = lfsv_s.template getChild<wetting>();
                const WettingPhase& w_space_n = lfsv_n.template getChild<wetting>();
                typedef typename WettingPhase::template Child<0>::Type WettingCellUnknowns;
                const WettingCellUnknowns& w_cell_space_s = w_space_s.template getChild<0>();
                const WettingCellUnknowns& w_cell_space_n = w_space_n.template getChild<0>();
                typedef typename WettingPhase::template Child<1>::Type WettingFaceUnknowns;
                const WettingFaceUnknowns& w_face_space_s = w_space_s.template getChild<1>();

                // local centers of adjacent cells
                GeometryType gt_s = ig.inside()->geometry().type();
                FieldVector<ctype,dim> localcenter_s = GenericReferenceElements<ctype,dim>::general(gt_s).position(0,0);
                GeometryType gt_n = ig.outside()->geometry().type();
                FieldVector<ctype,dim> localcenter_n = GenericReferenceElements<ctype,dim>::general(gt_n).position(0,0);

                // extract wetting phase saturation and non-wetting phase pressure
                rtype S_w_s = x_s(w_cell_space_s,0);
                rtype S_w_n = x_n(w_cell_space_n,0);
                rtype p_n_s = x_s(n_cell_space_s,0);
                rtype p_n_n = x_n(n_cell_space_n,0);

                // get capillary pressure
                rtype p_c_s = data.pc(*ig.inside(), localcenter_s, S_w_s);
                rtype p_c_n = data.pc(*ig.outside(), localcenter_n, S_w_n);

                // calculate fluxes (not including mobilities so far)
                unsigned e = ig.intersectionIndex();
                rtype u_n = 0.0;
                rtype u_w = 0.0;
                for(unsigned int f = 0, m = e*cell.num_faces; f < cell.num_faces; ++f, ++m)
                {
		  u_n += W[m] * (p_n_s - x_s(n_face_space_s,f));
		  u_w += W[m] * (p_n_s - p_c_s - x_s(w_face_space_s,f));
                }

                // gravity term
                rtype grav_flux = gravity_flux(e);
                u_n += data.rho_g(*ig.inside(), localcenter_s, p_n_s) * grav_flux;
                u_w += data.rho_l(*ig.inside(), localcenter_s, p_n_s - p_c_s) * grav_flux;

                // flux continuity on the faces
                r_s.accumulate(n_face_space_s,e,-u_n);
		r_s.accumulate(w_face_space_s,e,-u_w);

                if (init_mode)
                    return;

                // determine upwind mobility of non-wetting phase
                rtype S_w_s_upw, S_w_n_upw, S_w_s_downw, S_w_n_downw;
                if (u_n > 0.0)
                {
                    S_w_s_upw = S_w_s;
                    //S_w_n_upw = data.s_l(*ig.outside(), localcenter_n, p_c_s);
                    S_w_n_upw = S_w_s;
                    S_w_s_downw = S_w_n;
                    S_w_n_downw = S_w_n;
                }
                else
                {
                    //S_w_s_upw = data.s_l(*ig.inside(), localcenter_s, p_c_n);
                    S_w_s_upw = S_w_n;
                    S_w_n_upw = S_w_n;
                    S_w_s_downw = S_w_s;
                    S_w_n_downw = S_w_s;
                }
                rtype lambda_n_s = data.kr_g(*ig.inside(), localcenter_s, 1.0 - S_w_s_upw)
                    / data.mu_g(*ig.inside(), localcenter_s, p_n_s);
                rtype lambda_n_n = data.kr_g(*ig.outside(), localcenter_n, 1.0 - S_w_n_upw)
                    / data.mu_g(*ig.outside(), localcenter_n, p_n_n);
                rtype lambda_n = std::sqrt(lambda_n_s * lambda_n_n);
                // rtype ratio = std::abs((1.0 - S_w_s_upw) / (1.0 - S_w_s_downw));
                // if (ratio < 2.1)
                // {
                //     rtype lambda_avg = data.kr_g(*ig.inside(), localcenter_s, 1.0 - 0.5*(S_w_s_downw+S_w_s_upw))
                //         / data.mu_g(*ig.inside(), localcenter_s, p_n_s);
                //     if (ratio > 1.1)
                //         lambda_n = 1.0 * ((ratio-1.1) * lambda_n +
                //                           (2.1-ratio) * lambda_avg);
                //     else
                //         lambda_n = lambda_avg;
                // }

                // determine upwind mobility of wetting phase
                if ((u_w > 0.0) != (u_n > 0.0))   // only necessary if not the same as before
                {
                    if (u_w > 0.0)
                    {
                        S_w_s_upw = S_w_s;
                        S_w_n_upw = S_w_s;
                        //S_w_n_upw = data.s_l(*ig.outside(), localcenter_n, p_c_s);
                        S_w_s_downw = S_w_n;
                        S_w_n_downw = S_w_n;
                    }
                    else
                    {
                        S_w_s_upw = S_w_n;
                        //S_w_s_upw = data.s_l(*ig.inside(), localcenter_s, p_c_n);
                        S_w_n_upw = S_w_n;
                        S_w_s_downw = S_w_s;
                        S_w_n_downw = S_w_s;
                    }
                }
                // ratio = std::abs(S_w_s_upw / S_w_s_downw);
                rtype lambda_w_s = data.kr_l(*ig.inside(), localcenter_s, S_w_s_upw)
                    / data.mu_l(*ig.inside(), localcenter_s, p_n_s - p_c_s);
                rtype lambda_w_n = data.kr_l(*ig.outside(), localcenter_n, S_w_n_upw)
                    / data.mu_l(*ig.outside(), localcenter_n, p_n_n - p_c_n);
                rtype lambda_w = std::sqrt(lambda_w_s * lambda_w_n);
                // if (ratio < 2.1)
                // {
                //     rtype lambda_avg = data.kr_l(*ig.inside(), localcenter_s, 0.5*(S_w_s_downw+S_w_s_upw))
                //         / data.mu_l(*ig.inside(), localcenter_s, p_n_s - p_c_s);
                //     if (ratio > 1.1)
                //         lambda_w = 1.0 * ((ratio-1.1) * lambda_w +
                //                           (2.1-ratio) * lambda_avg);
                //     else
                //         lambda_w = lambda_avg;
                // }

                // apply fluxes to cell balance equation
                r_s.accumulate(n_cell_space_s,0,timestep * lambda_n * u_n);
		r_s.accumulate(w_cell_space_s,0,timestep * lambda_w * u_w);
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                const LFSV& lfsv_s, R& r_s) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
                const NonWettingPhase& n_space_s = lfsv_s.template getChild<nonwetting>();
                typedef typename NonWettingPhase::template Child<0>::Type NonWettingCellUnknowns;
                const NonWettingCellUnknowns& n_cell_space_s = n_space_s.template getChild<0>();
                typedef typename NonWettingPhase::template Child<1>::Type NonWettingFaceUnknowns;
                const NonWettingFaceUnknowns& n_face_space_s = n_space_s.template getChild<1>();

                typedef typename LFSV::template Child<wetting>::Type WettingPhase;
                const WettingPhase& w_space_s = lfsv_s.template getChild<wetting>();
                typedef typename WettingPhase::template Child<0>::Type WettingCellUnknowns;
                const WettingCellUnknowns& w_cell_space_s = w_space_s.template getChild<0>();
                typedef typename WettingPhase::template Child<1>::Type WettingFaceUnknowns;
                const WettingFaceUnknowns& w_face_space_s = w_space_s.template getChild<1>();

                // local center of intersection
                unsigned e = ig.intersectionIndex();
                GeometryType gt = ig.intersection().type();
                FieldVector<ctype,dim-1> face_center = GenericReferenceElements<ctype,dim-1>::general(gt).position(0,0);

                // local center of adjacent cell
                GeometryType gt_s = ig.inside()->geometry().type();
                FieldVector<ctype,dim> localcenter_s = GenericReferenceElements<ctype,dim>::general(gt_s).position(0,0);

                // extract wetting phase saturation and non-wetting phase pressure
                rtype S_w_s = x_s(w_cell_space_s,0);
                rtype p_n_s = x_s(n_cell_space_s,0);

                // get capillary pressure
                rtype p_c_s = data.pc(*ig.inside(), localcenter_s, S_w_s);

                // calculate fluxes (not including mobilities so far)
                rtype u_n = 0.0;
                rtype u_w = 0.0;
                for(unsigned int f = 0, m = e*cell.num_faces; f < cell.num_faces; ++f, ++m)
                {
		  u_n += W[m] * (p_n_s - x_s(n_face_space_s,f));
		  u_w += W[m] * (p_n_s - p_c_s - x_s(w_face_space_s,f));
                }

                // gravity term
                rtype grav_flux = gravity_flux(e);
                u_n += data.rho_g(*ig.inside(), localcenter_s, p_n_s) * grav_flux;
                u_w += data.rho_l(*ig.inside(), localcenter_s, p_n_s - p_c_s) * grav_flux;

                // consider boundary conditions
                if (data.bc_g(ig.intersection(), face_center, time) == 0) // Neumann boundary
                {
		  r_s.accumulate(n_face_space_s,e,-u_n);

                    if (!init_mode)
                        // apply flux to cell balance equation
		      r_s.accumulate(n_cell_space_s,0,timestep * u_n);
                }
                else // Dirichlet boundary
                {
		  r_s.accumulate(n_face_space_s,e,1e-10*x_s(n_face_space_s,e));

                    if (!init_mode)
                    {
                        // determine mobility
                        rtype lambda_n = data.kr_g(*ig.inside(), localcenter_s, 1.0 - S_w_s)
                            / data.mu_g(*ig.inside(), localcenter_s, p_n_s);
                        // apply flux to cell balance equation
                        r_s.accumulate(n_cell_space_s,0,timestep * lambda_n * u_n);
                    }
                }
                if (data.bc_l(ig.intersection(), face_center, time) == 0) // Neumann boundary
                {
		  r_s.accumulate(w_face_space_s,e,-u_w);

                    if (!init_mode)
                        // apply flux to cell balance equation
		      r_s.accumulate(w_cell_space_s,0,timestep * u_w);
                }
                else // Dirichlet boundary
                {
		  r_s.accumulate(w_face_space_s,e,1e-10*x_s(w_face_space_s,e));
		  // determine mobility
                    if (!init_mode)
                    {
                        rtype lambda_w = data.kr_l(*ig.inside(), localcenter_s, S_w_s)
                            / data.mu_l(*ig.inside(), localcenter_s, p_n_s - p_c_s);
                        // apply flux to cell balance equation
                        r_s.accumulate(w_cell_space_s,0,timestep * lambda_w * u_w);
                    }
                }
            }

            // lambda **********************************************

            template<typename EG, typename LFSV, typename R>
            void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
                const NonWettingPhase& n_space = lfsv.template getChild<nonwetting>();
                typedef typename NonWettingPhase::template Child<0>::Type NonWettingCellUnknowns;
                const NonWettingCellUnknowns& n_cell_space = n_space.template getChild<0>();
                typedef typename LFSV::template Child<wetting>::Type WettingPhase;
                const WettingPhase& w_space = lfsv.template getChild<wetting>();
                typedef typename WettingPhase::template Child<0>::Type WettingCellUnknowns;
                const WettingCellUnknowns& w_cell_space = w_space.template getChild<0>();

                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = GenericReferenceElements<ctype,dim>::general(gt).position(0,0);

		typedef typename LFSV::Traits::GridFunctionSpaceType::Traits::BackendType B;
                rtype S_w_old = B::access(x_old, w_cell_space.globalIndex(0));

                if (!init_mode)
                {
		  r.accumulate(n_cell_space,0,-cell.volume *
			       (phi * (1.0 - S_w_old) + data.q_g(eg.entity(), localcenter, time)));
		  r.accumulate(w_cell_space,0,-cell.volume *
		      (phi * S_w_old + data.q_l(eg.entity(), localcenter, time)));
                }
                else
                {
                    rtype p_n_old = B::access(x_old, n_cell_space.globalIndex(0));
                    r.accumulate(n_cell_space,0,-p_n_old);
		    r.accumulate(w_cell_space,0,-S_w_old);
                }
            }

            template<typename IG, typename LFSV, typename R>
            void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<nonwetting>::Type NonWettingPhase;
                const NonWettingPhase& n_space = lfsv.template getChild<nonwetting>();
                typedef typename NonWettingPhase::template Child<1>::Type NonWettingFaceUnknowns;
                const NonWettingFaceUnknowns& n_face_space = n_space.template getChild<1>();
                typedef typename LFSV::template Child<wetting>::Type WettingPhase;
                const WettingPhase& w_space = lfsv.template getChild<wetting>();
                typedef typename WettingPhase::template Child<1>::Type WettingFaceUnknowns;
                const WettingFaceUnknowns& w_face_space = w_space.template getChild<1>();

                // local index of current face
                unsigned int e = ig.intersectionIndex();

                GeometryType gt = ig.intersection().type();
                FieldVector<ctype,dim-1> center = GenericReferenceElements<ctype,dim-1>::general(gt).position(0,0);

                if (data.bc_g(ig.intersection(), center, time) == 0) // Neumann boundary
		  r.accumulate(n_face_space,e,cell.face_areas[e] * data.j_g(ig.intersection(), center, time));
                else // Dirichlet boundary
		  r.accumulate(n_face_space,e,-1e-10*data.g_g(ig.intersection(), center, time));

                if (data.bc_l(ig.intersection(), center, time) == 0) // Neumann boundary
		  r.accumulate(w_face_space,e,cell.face_areas[e] * data.j_l(ig.intersection(), center, time));
                else // Dirichlet boundary
		  r.accumulate(w_face_space,e,-1e-10*data.g_l(ig.intersection(), center, time));
            }

        private:
            rtype gravity_flux(unsigned e) const
            {
                rtype grav_flux = 0.0;
                for (unsigned int i = 0, m = e*dim; i < dim; ++i, ++m)
                    grav_flux += K_gravity[i] * cell.N[m];
                grav_flux *= cell.face_areas[e];
                return grav_flux;
            }

            const Data& data;
            const V& x_old;
            const WBuilder wbuilder;
            rtype time, timestep;
            bool init_mode;
            mutable unsigned cell_cache_id;
            mutable MimeticCellProperties<ctype,dim> cell;
            mutable typename Data::Traits::RangeType K_gravity;
            mutable std::vector<rtype> W;
            mutable rtype phi;
        };
    }
}

#endif
