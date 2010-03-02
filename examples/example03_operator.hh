#include"example02_operator.hh"
#include<dune/pdelab/localoperator/idefault.hh>

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *   \nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam B a function indicating the type of boundary condition
 */
template<class B>
class Example03LocalOperator :
  public Example02LocalOperator<B>,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  Example03LocalOperator (const B& b_, unsigned int intorder_=2)
    : Example02LocalOperator<B>(b_,intorder_), time(0.0)
  {}

  //! set time for subsequent evaluation
  void setTime (double t) {time = t;}

private:
  double time;
};
