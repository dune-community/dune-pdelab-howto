// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Abstract Parameter Interface required by the ParameterFactory class to create parameter classes on the fly
*/
#ifndef DUNE_PARAMETER_BASE_HH
#define DUNE_PARAMETER_BASE_HH

template<typename GV, typename RF>
class ParameterBase
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  virtual std::string name() const = 0;

  //! tensor diffusion coefficient
  virtual typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const = 0;

  //! velocity field
  virtual typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const = 0;

  //! sink term
  virtual typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const = 0;

  //! source term
  virtual typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const = 0;

  //! boundary condition type function
  virtual BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const = 0;

  //! Dirichlet boundary condition value
  virtual typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const = 0;

  //! Neumann boundary condition
  virtual typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const = 0;

  //! outflow boundary condition
  virtual typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const = 0;

};

#endif // DUNE_PARAMETER_BASE_HH
