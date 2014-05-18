// ParallelUnionFindImpl.h - declaration of the ParallelUnionFindImpl class which is a part of the Bridge pattern.
// Through the Template Method pattern defines the general sceleton
// of the parallel union-find algorithm described in:
// C.Harrison et al. Eurographics Symposium on Parallel Graphics and Visualization, 2011.

#ifndef PARALLEL_UNION_FIND_IMPL
#define PARALLEL_UNION_FIND_IMPL

#include <mpi.h>
#include <string>
#include "WeightedUnionFind.h"

struct DecompositionInfo
{
    DecompositionInfo() : numOfProc(0), myRank(-1), domainWidth(1), domainHeight(1) {}

    //DecompositionInfo(const DecompositionInfo& di): 
    // numOfProc(di.numOfProc), myRank(di.myRank), domainWidth(di.domainWidth), domainHeight(di.domainHeight) {}

    int numOfProc;
    int myRank;
    std::size_t domainWidth; // Width of the domain per processor, in grid points.
    std::size_t domainHeight;// Height of the domain per processor, in grid points.
};

class ParallelUnionFindImpl
{
public:
    ParallelUnionFindImpl(const DecompositionInfo& info);
    virtual ~ParallelUnionFindImpl(void);

    virtual void analyze(void);  // Template method pattern: specify a skeleton for the algorithm.

protected:
    // These are the fixed steps of the algorithm called from the analyze().
    virtual void runLocalUnionFind(void) = 0;
    virtual void constructGlobalLabeling(void) = 0;
    virtual void mergeLabelsAcrossProcessors(void) = 0;
    virtual void performFinalLabelingOfClusters(void) = 0;

protected:
    DecompositionInfo mDecompositionInfo;
};

#endif // PARALLEL_UNION_FIND_IMPL
