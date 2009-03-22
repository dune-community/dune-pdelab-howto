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
#include"dirichletboundaryconditionfunction.hh"
#include"laplacedirichletop.hh"

template<typename GV, typename FEM, typename CON> 
void laplacedirichlet (const GV& gv, const FEM& fem, 
					   int qorder, std::string filename)
{
  typedef typename GV::Grid::ctype D; // domain type                  /*@\label{lapdriver:FirstKnown}@*/
  typedef double R;                   // range type

  // set up grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
	Dune::PDELab::ISTLVectorBackend<1> > GFS; // Use ISTL Vector      /*@\label{lapdriver:ISTLBackend}@*/
  GFS gfs(gv,fem);                    // make grid function space

  typedef typename GFS::template ConstraintsContainer<R>::Type T;
  T t;                                // container for transformation
  B<GV> b(gv);                        // boundary condition function
  Dune::PDELab::constraints(b,gfs,t); // fill container

  typedef typename GFS::template VectorContainer<R>::Type X;
  X x(gfs,0.0);                       // make coefficient vector

  U<GV,R> u(gv);                      // analytic function object
  Dune::PDELab::interpolate(u,gfs,x); // interpolate x from u
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x); // clear interior  /*@\label{lapdriver:LastKnown}@*/

  // set up operator
  LaplaceDirichlet lop(qorder);       // local operator              /*@\label{lapdriver:ResEval}@*/
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,
    LaplaceDirichlet,T,T,
	Dune::PDELab::ISTLBCRSMatrixBackend<1,1> > GOS;                  /*@\label{lapdriver:GOS}@*/
  GOS gos(gfs,t,gfs,t,lop);           // global operator 

  // solve problem with matrix
  typedef typename GOS::template MatrixContainer<R>::Type M;         /*@\label{lapdriver:MatrixType}@*/
  M m(gos); m = 0.0;                  // matrix representation       /*@\label{lapdriver:MatrixSetup}@*/
  gos.jacobian(x,m);                  // assemble stiffness matrix   /*@\label{lapdriver:JacoAssemble}@*/

  Dune::MatrixAdapter<M,X,X> op(m);   // ISTL operator from matrix   /*@\label{lapdriver:ISTLFirst}@*/
  Dune::SeqSSOR<M,X,X> ssor(m,1,1.0); // a preconditioner
  Dune::CGSolver<X> solver(op,ssor,1E-10,5000,2); // CG solver
  Dune::InverseOperatorResult stat;   // status report of solver     /*@\label{lapdriver:ISTLLast}@*/

  X r(gfs,0.0);                       // residual vector             
  gos.residual(x,r);                  // compute residual            /*@\label{lapdriver:Residual}@*/
  X z(gfs,0.0);                       // update vector
  solver.apply(z,r,stat);             // solve for update            /*@\label{lapdriver:UpdateEq}@*/
  x -= z;                             // apply update                /*@\label{lapdriver:Update}@*/

  // solve problem without matrix ("on the fly")
  Dune::PDELab::OnTheFlyOperator<X,X,GOS> opb(gos);// operator apply /*@\label{lapdriver:OnTheFlyOp}@*/
  Dune::Richardson<X,X> richardson(1.0);           // no precond.    /*@\label{lapdriver:Richardson}@*/
  Dune::CGSolver<X> solverb(opb,richardson,1E-10,5000,2); // CG      /*@\label{lapdriver:NewCG}@*/
  
  Dune::PDELab::interpolate(u,gfs,x); // reinitialize x              /*@\label{lapdriver:AltSolverFirst}@*/
  Dune::PDELab::set_nonconstrained_dofs(t,0.0,x);  // clear interior
  r = 0.0;                            // clear residual
  gos.residual(x,r);                  // compute residual
  z = 0.0;                            // clear update
  solverb.apply(z,r,stat);            // solve for update 
  x -= z;                             // apply update                /*@\label{lapdriver:AltSolverLast}@*/

  // show result
  typedef Dune::PDELab::DiscreteGridFunction<GFS,X> DGF;             /*@\label{lapdriver:VTKFirst}@*/
  DGF dgf(gfs,x);                     // make a grid function
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0); // plot result
  vtkwriter.addVertexData(new Dune::PDELab::
						  VTKGridFunctionAdapter<DGF>(dgf,"u"));
  vtkwriter.write(filename,Dune::VTKOptions::ascii);                 /*@\label{lapdriver:VTKLast}@*/
}
