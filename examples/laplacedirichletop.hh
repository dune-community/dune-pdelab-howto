#include<vector>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

class LaplaceDirichlet
: public Dune::PDELab::NumericalJacobianApplyVolume<LaplaceDirichlet>, /*@\label{laplace:JacApply}@*/
  public Dune::PDELab::NumericalJacobianVolume<LaplaceDirichlet>,      /*@\label{laplace:Jac}@*/
  public Dune::PDELab::FullVolumePattern,                              /*@\label{laplace:Pattern}@*/
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };    /*@\label{laplace:FirstFlag}@*/

  // residual assembly flags
  enum { doAlphaVolume = true };     /*@\label{laplace:LastFlag}@*/

  LaplaceDirichlet (int qorder_) : qorder(qorder_) {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, 
		   typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, /*@\label{laplace:AlphaVolume}@*/
					 const LFSV& lfsv, R& r) const
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
	const Dune::QuadratureRule<DF,dim>& 
	  rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

	// loop over quadrature points
	for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
		   it=rule.begin(); it!=rule.end(); ++it)
	  {
		// evaluate gradient of shape functions (assume Galerkin)
		std::vector<JacobianType> js(lfsu.size());
		lfsu.localFiniteElement().localBasis().evaluateJacobian(
		  it->position(),js);

		// transform gradient to real element
		const Dune::FieldMatrix<DF,dimw,dim> jac = 
		  eg.geometry().jacobianInverseTransposed(it->position());
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
		RF factor = it->weight() * eg.geometry().integrationElement(
		                                          it->position());
		for (int i=0; i<lfsu.size(); i++)
		  r[i] += (gradu*gradphi[i])*factor;
	  }
  }
private:
  int qorder;
};

