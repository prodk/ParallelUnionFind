ParallelUnionFind
_________________________________________________________
C++ implementation of the parallel union-find algorithm for distributed-memory machines 
(can be used on multi-core CPUs as well as multi-node supercomputers).

_________________________________________________________
The algorithm is described in:
C.Harrison, H.Childs, K.P.Gaither,
Data-Parallel Mesh Connected Components Labeling and Analysis,
Eurographics Symposium on Parallel Graphics and Visualisation, 2011, (eds. T.Kuhlen, R.Pajarola, K.Zhou).

_________________________________________________________
(c) Mykola Prodanov, spring-autumn 2014, Odessa, Ukraine.

_________________________________________________________
How to use:
1) fill in the fields of the DecompositionInfo structure: 
provide the mesh size, the number of processors, a pointer to the data etc.;
2) create a ParallelUnionFind object of the necessary configuration using DecompositionInfo structure;
3) call analyze() of the ParallelUnionFind object;
4) print the necessary output via the corresponding interface of the ParallelUnionFind object.

_________________________________________________________
Restrictions/assumptions of the current implementation:

1) only a 2D domain decomposition in the form of parallel stripes is supported ("2DStripes" configuration);
2) message passing using MPI is assumed;
3) no binary space partitioning (BSP) for optimizing load-balancing is used
(we assume that the mesh is uniformly distributed over the space);
4) we assume that the data is already present on every node involved in the calculations.
You have to provide a pointer to it.
5) the data is assumed to be a 1D array. You have to convert your data into this format.

_________________________________________________________
To extend the implementation, e.g. for supporting a 3D domain decomposition,
do the following:

1) make your own implementation class via inheriting 
from the Template Method pattern abstract base class ParallelUnionFindImpl;
2) introduce your configuration, e.g. "3DCubeDecomposition";
3) add to the factory constructing the object of your configuration;
4) override basic steps defined in ParallelUnionFindImpl, use the helper serial WeightedUnionFind class.