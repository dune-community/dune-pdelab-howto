#ifndef CG_STOKES_INITIAL_HH
#define CG_STOKES_INITIAL_HH

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

// constraints parameter class for selecting boundary condition type
class BCTypeParam_HagenPoiseuille
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParam_HagenPoiseuille() {}

  template<typename I>
  inline void evaluate (
    const I & intersection,   /*@\label{bcp:name}@*/
    const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
    BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
        xg = intersection.geometry().global( coord );
    if( xg[0] < 1e-6 )
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
};



// constraints parameter class for selecting boundary condition type
class BCTypeParam_TU
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParam_TU() {}

  template<typename I>
  inline void evaluate (
    const I & intersection,   /*@\label{bcp:name}@*/
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




// constraints parameter class for selecting boundary condition type
class BCTypeParam_LU
{
public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  BCTypeParam_LU() {}

  template<typename I>
  inline void evaluate (
    const I & intersection,   /*@\label{bcp:name}@*/
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





// constraints parameter class for selecting boundary condition type
template<typename GV>
class BCTypeParam_PressureDrop
{
private:
  typedef typename GV::ctype DFT;
  const DFT length;
  const DFT origin;
  const int direction;

public:
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };


  BCTypeParam_PressureDrop (const DFT l_, const DFT o_, const int d_)
    : length(l_), origin(o_), direction(d_)
  {
  }

  template<typename I>
  inline void evaluate (
    const I & intersection,   /*@\label{bcp:name}@*/
    const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
    BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
      xg = intersection.geometry().global( coord );

    if(xg[direction]-origin < 1e-6 || xg[direction]-origin > length-1e-6)
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }

};

template<typename GV, typename RF, int dim>
class HagenPoiseuilleVelocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  HagenPoiseuilleVelocity<GV,RF,dim> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HagenPoiseuilleVelocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  HagenPoiseuilleVelocity(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    RF r = 0;

    for(int i=1; i<dim; ++i){
      r += (x[i]-0.5)*(x[i]-0.5);
      y[i] = 0;
    }
    r = sqrt(r);
    y[0] = 0.25 - r*r;
  }

};


template<typename GV, typename RF, int dim>
class TU_Velocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  TU_Velocity<GV,RF,dim> >
{
  RF time;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,TU_Velocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  TU_Velocity(const GV & gv) : BaseT(gv) {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    RF r = 0;

    // 2d:
    //y[1] = 0;
    //r -= (x[1]-0.5)*(x[1]+0.5);

    for(int i=0;i<dim;i++){
      y[i] = 0;
      if(i>0){
        r -= (x[i]-0.5)*(x[i]+0.5);
      }
    }

    if(x[0] < 1e-6){
      y[0] = r*(1.0 - std::cos(2.0*3.1415*time));
    }
    else{
      y[0]=0;
    }
  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};



template<typename GV, typename RF, int dim>
class LU_Velocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
  LU_Velocity<GV,RF,dim> >
{
private:
  RF time;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,LU_Velocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  LU_Velocity(const GV & gv) : BaseT(gv) {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    RF r = 0;

    y[1] = 0;
    r += x[1]*(1.0-x[1]);

    if(x[0] < -1.0+1e-6){
      y[0] = r*(1.0 - std::cos(2.0*3.1415*time));
    }
    else{
      y[0]=0;
    }
  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};




template<typename GV, typename RF>
class HagenPoiseuilleVelocity<GV,RF,3> :
  public Dune::PDELab::AnalyticGridFunctionBase<
  Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,3>,
  HagenPoiseuilleVelocity<GV,RF,3> >
{
public:
  enum {dim = 3};
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, HagenPoiseuilleVelocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  HagenPoiseuilleVelocity(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {

    RF r = 0;
    for(int i=1; i<dim; ++i){
      r += (x[i])*(x[i]);
      y[i] = 0;
    }
    r = sqrt(r);
    y[0] = 0.25 - r*r;

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
