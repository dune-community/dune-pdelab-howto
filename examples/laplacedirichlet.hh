#include<vector>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>

class LaplaceDirichlet 
  : public Dune::PDELab::NumericalJacobianApplyVolume<LaplaceDirichlet>,
	public Dune::PDELab::NumericalJacobianVolume<LaplaceDirichlet>,
	public Dune::PDELab::FullVolumePattern
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = false };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaSkeleton = false };
  enum { doAlphaBoundary = false };
  enum { doLambdaVolume = false };
  enum { doLambdaSkeleton = false };
  enum { doLambdaBoundary = false };

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
	  Traits::LocalBasisType::Traits::JacobianType JacobianType;

	// dimensions
	const int dim = EG::Geometry::dimension;
	const int dimw = EG::Geometry::dimensionworld;

	// select quadrature rule
	Dune::GeometryType gt = eg.geometry().type();
	const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,3);

	// loop over quadrature points
	for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
	  {
		// evaluate gradient of shape functions (we assume Galerkin method lfsu=lfsv)
		std::vector<JacobianType> js(lfsu.size());
		lfsu.localFiniteElement().localBasis().evaluateJacobian(it->position(),js);

		// transform gradient to real element
		const Dune::FieldMatrix<DF,dimw,dim> jac = eg.geometry().jacobianInverseTransposed(it->position());
		std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
		for (int i=0; i<lfsu.size(); i++)
		  {
			gradphi[i] = 0.0;
			jac.umv(js[i][0],gradphi[i]);
		  }

		// compute gradient of u
		Dune::FieldVector<RF,dim> gradu(0.0);
		for (int i=0; i<lfsu.size(); i++)
		  gradu.axpy(x[i],gradphi[i]);

		// integrate grad u * grad phi_i
		RF factor = it->weight() * eg.geometry().integrationElement(it->position());
		for (int i=0; i<lfsu.size(); i++)
		  r[i] += (gradu*gradphi[i])*factor;
	  }
  }
};

