#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include"q1dirichletboundaryconditionfunction.hh"
#include"laplacedirichlet.hh"

template<typename GV>
void q1laplacedirichlet (const GV& gv)
{
  typedef typename GV::Grid::ctype D; // domain type
  typedef double R;                   // range type

  // set up function space
  Q1LocalFiniteElementMap<D,R> fem;   // maps entity to finite element

  typedef Dune::PDELab::GridFunctionSpace<GV,
	Q1LocalFiniteElementMap<D,R>,Q1Constraints,
	Dune::PDELab::ISTLVectorBackend<1> > GFS; // Use ISTL Vector    
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template ConstraintsContainer<R>::Type T;
  T t;                                // container for transformation
  B<GV> b(gv);                        // boundary condition function
  Dune::PDELab::constraints(b,gfs,t); // fill container

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // make analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // interpolate x from u
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x); // clear interior

  // set up operator
  LaplaceDirichlet lop;               // local operator
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LaplaceDirichlet,T,T,
	Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;
  GOS gos(gfs,t,gfs,t,lop);           // global operator

  // solve problem with matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;
  M m(gos);                           // matrix representation
  m = 0.0;                            // clear matrix
  gos.jacobian(x,m);                  // assemble stiffness matrix

  Dune::MatrixAdapter<M,X,X> op(m);   // ISTL operator from matrix
  Dune::SeqSSOR<M,X,X> ssor(m,1,1.0); // a preconditioner
  Dune::CGSolver<X> solver(op,ssor,1E-10,5000,2); // CG solver
  Dune::InverseOperatorResult stat;   // status report of solver

  X r(gfs,0.0);                       // residual vector
  gos.residual(x,r);                  // compute residual
  X z(gfs,0.0);                       // update vector
  solver.apply(z,r,stat);             // solve for update 
  x -= z;                             // apply update

  // solve problem without matrix ("on the fly")
  Dune::PDELab::OnTheFlyOperator<X,X,GOS> opb(gos);// operator apply
  Dune::Richardson<X,X> richardson(1.0);           // no precond.
  Dune::CGSolver<X> solverb(opb,richardson,1E-10,5000,2); // CG
  
  Dune::PDELab::interpolate(u,gfs,x); // reinitialize x
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x);  // clear interior
  r = 0.0;                            // clear residual
  gos.residual(x,r);                  // compute residual
  z = 0.0;                            // clear update
  solverb.apply(z,r,stat);            // solve for update 
  x -= z;                             // apply update

  // show result
  typedef Dune::PDELab::DiscreteGridFunction<GFS,X> DGF;
  DGF dgf(gfs,x);                     // make a grid function

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1); // plot result
  vtkwriter.addVertexData(new Dune::PDELab::
						  VTKGridFunctionAdapter<DGF>(dgf,"q1"));
  vtkwriter.write("q1laplace",Dune::VTKOptions::ascii); 
}
