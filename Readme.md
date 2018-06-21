Readme
======

MPI+Pthreads Speckle simulation. Master thesis project.

This is a parallel speckle simulation for a Condensed Matter Physics experiment.

The system solves an eigenvalues problem in parallel with a standard API to face
mkl, plasma and magma. The ensemble averaging parallelism uses a master-slave
approach with an extra thread in the master node to improve the resources usage.

The master slave implementation uses a common API that allows the code to be
reused. That code will be a different project in a while.

**Author: Jimmy Aguilar Mena**
