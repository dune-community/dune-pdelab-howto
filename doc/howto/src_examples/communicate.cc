  // get communication object from gridview
  auto comm = gridview.comm();
  // get rank from communication object
  int myRank = comm.rank();
  // get number of processes
  int numProcs = comm.size();
  // calculate global sum (using MPI_Reduce)
  double globalsum=comm.sum(localResult);
  // calculate global maximum (using MPI_Reduce)
  double globalmax=comm.max(localResult);
  // broadcast result
  comm.broadcast(&globalMax,1,0);
  // wait until all processes arrived here
  comm.barrier();
