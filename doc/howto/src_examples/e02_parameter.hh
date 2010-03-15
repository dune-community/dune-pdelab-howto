#ifndef __E02_PARAMETER_HH_
#define __E02_PARAMETER_HH_

// layout for codim0 data
template <int dim>
struct P0Layout
{
  bool contains(Dune::GeometryType gt)
  {
    if ( gt.dim() == dim ) return true;
    return false;
  }
};

//===============================================================
// Leiffel parameter
// The L is a building block of the Eiffel tower
//===============================================================

/**
 * \brief A function defining the diffusion coefficient
 */
template<typename GV, typename RF, typename PGMap>
class LeiffelDiffusion
  : public Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
        LeiffelDiffusion<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,
      Dune::FieldVector<RF,1> >, LeiffelDiffusion<GV,RF,PGMap> > BaseT;

  // constructor
  LeiffelDiffusion(const GV& gv_, const PGMap& pg_) : gv(gv_), mapper(gv), pg(pg_) {}

  // evaluate scalar diffusion parameter
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // retrieve element index and corresponding material parameter index
    const int ei              = mapper.map(e);
    const int physgroup_index = pg[ei];

    // evaluate physical group map and set values accordingly
    switch ( physgroup_index )
    {
      case 1  : y = 2.0; break;
      case 2  : y = 2.0;  break;
      default : y = 1.0;  break;
    }
  }

  inline const typename Traits::GridViewType& getGridView() { return gv; }

private:

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
};

/** \brief boundary grid function selecting boundary conditions
 * 0 means Neumann
 * 1 means Dirichlet
 */
template<typename GV, typename PGMap>
class LeiffelBCType : public Dune::PDELab::BoundaryGridFunctionBase<
        Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
        Dune::FieldVector<int,1> >,LeiffelBCType<GV,PGMap> >
{

public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view
  LeiffelBCType (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    int physgroup_index = pg[i.boundarySegmentIndex()];
    switch ( physgroup_index )
    {
      case 3  : y = 0; break; // Neumann
      case 4  : y = 1; break; // Dirichlet
      default : y = 0; break; // Neumann is default
    }
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

private:

  const GV&    gv;
  const PGMap& pg;
};

/**
 * \brief A function that defines Dirichlet boundary conditions
 *        AND its extension to the interior
 */
template<typename GV, typename RF, typename PGMap>
class LeiffelBCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           LeiffelBCExtension<GV,RF,PGMap> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  //! construct from grid view
  LeiffelBCExtension(const GV& gv_, const PGMap& pg_) : gv(gv_), mapper(gv), pg(pg_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    // retrieve element index and corresponding material parameter index
    const int ei              = mapper.map(e);
    const int physgroup_index = pg[ei];

    // find Dirichlet boundaries and handle them
    y = 0.0;
    for (IntersectionIterator is=gv.ibegin(e); is!=gv.iend(e); ++is)
    {
      if ( is->boundary() )
      {
        if ( physgroup_index == 4 ) y = 1.0;
        break;
      }
    }
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
};

/**
 * \brief A function defining Flux boundary conditions
 */
template<typename GV, typename RF, typename PGMap>
class LeiffelFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, LeiffelFlux<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  // constructor
  LeiffelFlux(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  // evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i, const typename Traits::DomainType& xlocal,
                       typename Traits::RangeType& y) const
  {
    int physGroupInd = pg[i.boundarySegmentIndex()];
    switch ( physGroupInd )
    {
      case 3  : y = 0.0; break; // some influx
      default : y = 0.0; break; // Neumann-0 BC
    }
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
};

#endif
