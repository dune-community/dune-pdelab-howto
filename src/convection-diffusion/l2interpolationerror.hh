// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_L2INTERPOLARION_ERROR_HH
#define DUNE_PDELAB_L2INTERPOLARION_ERROR_HH

#include <cmath>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

/**
 * \brief Calculate error in the \f$L_2\f$ norm for a given
 * analytical function \a u and a result \a x bound to a grid
 * functions space \a gfs.
 *
 * \param u Analytical function
 * \param gfs GridFunctionSpace
 * \param x Solution bound to the GridFunctionSpace
 */
template<class U, class GFS, class X>
double l2interpolationerror (const U& u, const GFS& gfs, X& x,
                             int qOrder=1)
{
  // constants and types
  typedef typename GFS::Traits::GridViewType GV;
  const int dim = GV::dimension;
  typedef typename GV::Traits::template Codim<0>::Iterator
    ElementIterator;
  typedef typename GFS::Traits::FiniteElementType::
    Traits::LocalBasisType::Traits FETraits;
  typedef typename FETraits::DomainFieldType D;
  typedef typename FETraits::RangeFieldType R;
  typedef typename FETraits::RangeType RangeType;

  // make local function space
  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename X::template ConstLocalView<LFSCache> XView;

  LFS lfs(gfs);
  LFSCache lfsCache(lfs);
  XView xView(x);
  // local coefficients
  std::vector<R> xLocal(lfs.maxSize());
  // shape function values
  std::vector<RangeType> basisLocal(lfs.maxSize());

  // loop over grid view
  double sum = 0.0;
  for (ElementIterator eIt = gfs.gridView().template begin<0>();
       eIt!=gfs.gridView().template end<0>(); ++eIt)
  {
    // bind local function space to element
    lfs.bind(*eIt);
    lfsCache.update();
    xView.bind(lfsCache);
    xView.read(xLocal);
    xView.unbind();

    // integrate over element using a quadrature rule
    Dune::GeometryType gt = eIt->geometry().type();
    const Dune::QuadratureRule<D,dim>& rule =
      Dune::QuadratureRules<D,dim>::rule(gt, qOrder);

    for (typename Dune::QuadratureRule<D, dim>::const_iterator qIt = rule.begin();
          qIt != rule.end(); ++qIt)
    {
      // evaluate solution bound to grid function space at integration point
      RangeType uGfs(0.0);
      lfs.finiteElement().localBasis().evaluateFunction(qIt->position(), basisLocal);
      for (unsigned int i = 0; i < lfs.size(); ++i)
      {
        uGfs.axpy(xLocal[i], basisLocal[i]);
      }

      // evaluate the analytic function at integration point
      RangeType uAnalytic;
      u.evaluate(*eIt, qIt->position(), uAnalytic);

      // accumulate error
      uGfs -= uAnalytic;
      sum += uGfs.two_norm2() * qIt->weight()*
        eIt->geometry().integrationElement(qIt->position());
    }
  }
  return std::sqrt(sum);
}

#endif // DUNE_PDELAB_L2INTERPOLARION_ERROR_HH
