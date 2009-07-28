// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_POISSONHDIV_HH
#define DUNE_PDELAB_POISSONHDIV_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {

	// a local operator for solving the Poisson equation
	//           - \Delta u = f         in \Omega, 
    //                    u = g         on \partial\Omega_D
    //  -\nabla u \cdot \nu = v\cdot\nu on \partial\Omega_N
	// with H(div) conforming (mixed) finite elements
    // F : grid function type giving f
    // B : grid function type selecting boundary condition
    // G : grid function type giving g
    template<typename F, typename B, typename G, int qorder_v=2, int qorder_p=1>
	class PoissonHDivConforming : public NumericalJacobianApplyVolume<PoissonHDivConforming<F,B,G,qorder_v,qorder_p> >,
                                  public NumericalJacobianVolume<PoissonHDivConforming<F,B,G,qorder_v,qorder_p> >,
                                  public FullVolumePattern,
                                  public LocalOperatorDefaultFlags
	{
	public:
      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      PoissonHDivConforming (const F& f_, const B& b_, const G& g_)
        : f(f_), b(b_), g(g_)
      {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
        // select the two components
        typedef typename LFSU::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsu.template getChild<0>();
        typedef typename LFSU::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsu.template getChild<1>();

		// domain and range field type
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType VelocityJacobianType;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType VelocityRangeType;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        const int dimw = EG::Geometry::dimensionworld;

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = eg.geometry().integrationElement(pos);

        // \sigma\cdot v term
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& vrule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_v);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=vrule.begin(); it!=vrule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateFunction(it->position(),vbasis);

            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (int i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // compute sigma
            VelocityRangeType sigma; sigma = 0.0;
            for (int i=0; i<velocityspace.size(); i++)
              sigma.axpy(x[velocityspace.localIndex(i)],vtransformedbasis[i]);

            // integrate sigma * phi_i
            RF factor = it->weight() / det;
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] += (sigma*vtransformedbasis[i])*factor;
          }

        // u div v term and div sigma q term
        const Dune::QuadratureRule<DF,dim>& prule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=prule.begin(); it!=prule.end(); ++it)
          {
            // evaluate shape functions at ip (this is a Galerkin method)
            std::vector<VelocityJacobianType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateJacobian(it->position(),vbasis);
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.localFiniteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // compute u
            PressureRangeType u; u = 0.0;
            for (int i=0; i<pressurespace.size(); i++)
              u.axpy(x[pressurespace.localIndex(i)],pbasis[i]);

            // compute divergence of velocity basis functions on reference element
            std::vector<RF> divergence(velocityspace.size(),0.0);
            for (int i=0; i<velocityspace.size(); i++)
              for (int j=0; j<dim; j++) 
                divergence[i] += vbasis[i][j][j];

            // integrate sigma * phi_i
            RF factor = it->weight();
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] -= u*divergence[i]*factor;

            // compute divergence of sigma
            RF divergencesigma = 0.0;
            for (int i=0; i<velocityspace.size(); i++)
              divergencesigma += x[velocityspace.localIndex(i)]*divergence[i];

            // integrate div sigma * q
            for (int i=0; i<pressurespace.size(); i++)
              r[pressurespace.localIndex(i)] -= divergencesigma*pbasis[i]*factor;
          }
	  }

 	  // volume integral depending only on test functions
	  template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<1>::Type PressureSpace;
        const PressureSpace& pressurespace = lfsv.template getChild<1>();

		// domain and range field type
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename PressureSpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType PressureRangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder_p);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate shape functions 
            std::vector<PressureRangeType> pbasis(pressurespace.size());
            pressurespace.localFiniteElement().localBasis().evaluateFunction(it->position(),pbasis);

            // evaluate right hand side parameter function
            typename F::Traits::RangeType y;
            f.evaluate(eg.entity(),it->position(),y);

            // integrate f
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (int i=0; i<pressurespace.size(); i++)
              r[pressurespace.localIndex(i)] += y*pbasis[i]*factor;
          }
      }

      // boundary integral independen of ansatz functions
 	  template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type VelocitySpace;
        const VelocitySpace& velocityspace = lfsv.template getChild<0>();

		// domain and range field type
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename VelocitySpace::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType VelocityRangeType;

        // dimensions
        const int dim = IG::dimension;
        const int dimw = IG::dimensionworld;

        // evaluate transformation which must be linear
        Dune::FieldVector<DF,dim> pos; pos = 0.0;
        Dune::FieldMatrix<DF,dimw,dim> jac = ig.inside()->geometry().jacobianInverseTransposed(pos);
        jac.invert();
        RF det = ig.inside()->geometry().integrationElement(pos);

        // select quadrature rule
        Dune::GeometryType gtface = ig.geometryInInside().type();
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder_v);

        // loop over quadrature points and integrate normal flux
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
          {
            // evaluate boundary condition type
            typename B::Traits::RangeType bctype;
            b.evaluate(ig,it->position(),bctype);
 
            // skip rest if we are on Neumann boundary
            if (bctype<=0) continue;

            // position of quadrature point in local coordinates of element 
            Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

            // evaluate test shape functions 
            std::vector<VelocityRangeType> vbasis(velocityspace.size());
            velocityspace.localFiniteElement().localBasis().evaluateFunction(local,vbasis);
            
            // transform basis vectors
            std::vector<VelocityRangeType> vtransformedbasis(velocityspace.size());
            for (int i=0; i<velocityspace.size(); i++)
              {
                vtransformedbasis[i] = 0.0;
                jac.umtv(vbasis[i],vtransformedbasis[i]);
              }

            // evaluate Dirichlet boundary condition
            typename G::Traits::RangeType y;
            g.evaluate(*(ig.inside()),local,y);
            
            // integrate g v*normal
            RF factor = it->weight()*ig.geometry().integrationElement(it->position())/det;
            for (int i=0; i<velocityspace.size(); i++)
              r[velocityspace.localIndex(i)] += y*(vtransformedbasis[i]*ig.unitOuterNormal(it->position()))*factor;
          }
      }

    private:
      const F& f;
      const B& b;
      const G& g;
	};

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif
