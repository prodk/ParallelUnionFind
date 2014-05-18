// ParallelUnionFindImpl.cpp
// Implementation of the ParallelUnionFindImpl class.

#include "ParallelUnionFindImpl.h"

ParallelUnionFindImpl::ParallelUnionFindImpl(int numOfProc, int myRank):
mNumOfProc(numOfProc), mMyRank(myRank)
{
}

ParallelUnionFindImpl::~ParallelUnionFindImpl(void)
{
}

// Provides a skeleton for the parallel union find algorithm.
void ParallelUnionFindImpl::analyze()
{
    runLocalUnionFind();
    constructGlobalLabeling();
    mergeLabelsAcrossProcessors();
    performFinalLabelingOfClusters();
}
