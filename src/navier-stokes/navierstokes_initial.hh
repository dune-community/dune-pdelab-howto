#ifndef CG_STOKES_INITIAL_HH
#define CG_STOKES_INITIAL_HH

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

/** \brief Global Dirichlet boundary conditions

 */

// constraints parameter class for selecting boundary condition type
class BCTypeParamGlobalDirichlet
{
public :
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParamGlobalDirichlet() {}

  template<typename I>
  inline void evaluate (const I & intersection,   /*@\label{bcp:name}@*/
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    y = BC::VelocityDirichlet;
  }

};

/** \brief Boundary conditions for Hagen-Poiseuille flow.

 */

// constraints parameter class for selecting boundary condition type
class BCTypeParamHagenPoiseuille
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParamHagenPoiseuille() {}

  template<typename I>
  inline void evaluate (const I & intersection,   /*@\label{bcp:name}@*/
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    if( xg[0] > 1.0-1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
};

/** \brief Boundary conditions for two dimensional flow around a cylinder obstacle.

 * The settings come from the two dimensional Turek benchmark channel.
 * Compare with the Article: "Benchmark Computations of Laminar Flow Around a Cylinder"
 * by M. Schaefer and S. Turek.

 */

// constraints parameter class for selecting boundary condition type
class BCTypeParamTU
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParamTU() {}

  template<typename I>
  inline void evaluate (const I & intersection,   /*@\label{bcp:name}@*/
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    if( xg[0] > 2.2-1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
};

/** \brief Boundary conditions for three dimensional flow around a cylinder obstacle

 * The settings come from the three dimensional Turek benchmark channel.
 * Compare with the Article: "Benchmark Computations of 3D Laminar Flow Around a Cylinder with CFX, OpenFOAM and FEATFLOW"
 * by E. Bayraktar, O. Mierka and S. Turek.

 */

// constraints parameter class for selecting boundary condition type
class BCTypeParamTU3D
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParamTU3D() {}

  template<typename I>
  inline void evaluate (const I & intersection,   /*@\label{bcp:name}@*/
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    if( xg[0] > 2.5-1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
};

/** \brief Boundary conditions for Backward-Facing Step.

 */

// constraints parameter class for selecting boundary condition type
class BCTypeParamLU
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParamLU() {}

  template<typename I>
  inline void evaluate (const I & intersection,   /*@\label{bcp:name}@*/
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );
    if( xg[0] > 5.0-1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
};

/** \brief Dirichlet velocity function for Hagen-Poiseuille flow in a box.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.
    \tparam dim Dimension of the problem.

    * The Hagen-Poiseuille flow is solved on the dim-dimensional unitcube.
    * The maximum inflow velocity is normalized to 1.
    * In three dimensions, the two dimensional flow is trivially extended
    * to the z-direction.

    */

template<typename GV, typename RF, int dim>
class HagenPoiseuilleVelocityBox :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  HagenPoiseuilleVelocityBox<GV,RF,dim> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HagenPoiseuilleVelocityBox<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  HagenPoiseuilleVelocityBox(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y = 0;
    y[0] = x[1]*(1.0-x[1]);
    y[0] *= 4.0;
  }

};

/** \brief Dirichlet velocity function for the flow around a cylinder obstacle.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.
    \tparam dim The dimension of the problem.
    \param[in] meanflow Maximum inflow velocity.

    * This function is used to model the inflow.
    * It works both in two and three dimensions.

    */

template<typename GV, typename RF, int dim>
class TUVelocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  TUVelocity<GV,RF,dim> >
{
  RF time;
  const RF& meanflow_;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,TUVelocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  TUVelocity(const GV & gv, const RF& meanflow) :
    BaseT(gv)
    , meanflow_(meanflow)
  {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y = 0.0;
    RF r = 1.0;

    for(int i=1; i<dim; i++)
      r *= 4.0*x[i]*(0.41-x[i])/(0.41*0.41);

    if(x[0] < 1e-6) {
      y[0] = r*meanflow_;
      if(time <= 4.0)
        y[0] *= sin(M_PI*time/8.0);
    }
  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};

/** \brief Dirichlet velocity function for the Backward-Facing Step.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.
    \tparam dim The dimension of the problem.

    * The inflow velocity is normalized to 1.

    */

template<typename GV, typename RF, int dim>
class LUVelocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  LUVelocity<GV,RF,dim> >
{
private:
  RF time;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,LUVelocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  LUVelocity(const GV & gv) : BaseT(gv) {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    RF r = 0;

    y = 0.0;
    r += x[1]*(1.0-x[1]);

    if(x[0] < -1.0+1e-6) {
      y[0] = 4*r;
      if(time <= 2.0)
        y[0] *= sin(M_PI/4*time);
    }
  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};

/** \brief Dirichlet velocity function for three dimensional Hagen-Poiseuille flow.

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.

    * The Hagen-Poiseuille flow is solved in a pipe pointing in x-direction.
    * The cross section in the yz-plane is a circle with midpoint at the origin
    * and with radius 0.5.
    * The inflow velocity is normalized to 1.

    */

template<typename GV, typename RF>
class HagenPoiseuilleVelocityCylindrical :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
  HagenPoiseuilleVelocityCylindrical<GV,RF> >
{
public:
  enum {dim = 3};
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HagenPoiseuilleVelocityCylindrical<GV,RF> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  HagenPoiseuilleVelocityCylindrical(const GV & gv) : BaseT(gv)
  {
    assert(GV::dimension == dim);
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {

    RF r = 0;
    for(int i=1; i<dim; ++i){
      r += (x[i])*(x[i]);
      y[i] = 0;
    }
    r = sqrt(r);
    y[0] = 0.25 - r*r;
    y[0] *= 4;

    if(y[0]<0)
      y[0] = 0;
  }

};


template<typename GV, typename RF, std::size_t dim_range>
class ZeroVectorFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range>,
  ZeroVectorFunction<GV,RF,dim_range> >,
  public Dune::PDELab::InstationaryFunctionDefaults
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroVectorFunction> BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  ZeroVectorFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=0;
  }
};

template<typename GV, typename RF>
class ZeroScalarFunction
  : public ZeroVectorFunction<GV,RF,1>
{
public:

  ZeroScalarFunction(const GV & gv) : ZeroVectorFunction<GV,RF,1>(gv) {}

};

/** \brief Do-nothing boundary function for the Hagen-Poiseuille flow

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.

*/

// function for defining the flux boundary condition
template<typename GV, typename RF>
class HagenPoiseuilleZeroFlux
  : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  HagenPoiseuilleZeroFlux<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,HagenPoiseuilleZeroFlux<GV,RF> > BaseT;

  HagenPoiseuilleZeroFlux (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 0;
  }
};

/** \brief Do-nothing boundary function for the Backward-Facing Step and the flow around a cylinder

    \tparam GV The underlying grid view.
    \tparam RF The type of the range field.
    \param[in] p_ The pressure at the inflow.
    \param[in] l_ The length of the tube.
    \param[in] o_ The coordinate referring to the beginning of the tube.
    \param[in] d_ The canonical direction the tube is oriented.

*/

// function for defining the flux boundary condition
template<typename GV, typename RF>
class PressureDropFlux
  : public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
  PressureDropFlux<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,PressureDropFlux<GV,RF> > BaseT;

private:
  typedef typename Traits::DomainFieldType DFT;

  const DFT pressure;
  const DFT length;
  const DFT origin;
  const int direction;

public:
  PressureDropFlux (const GV& gv, const RF p_, const RF l_, const RF o_, const int d_)
    : BaseT(gv), pressure(p_), length(l_), origin(o_), direction(d_)
  {
#ifndef NDEBUG
    const int dim = GV::dimension;
    assert(direction >=0 && direction <dim);
#endif
  }

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    if(x[direction]-origin < 1e-6)
      y = pressure;
    else if(x[direction]-origin > length-1e-6)
      y = 0.;
  }
};


#endif
