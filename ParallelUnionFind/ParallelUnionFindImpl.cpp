// ParallelUnionFindImpl.cpp
// Implementation of the ParallelUnionFindImpl class.

#include "ParallelUnionFindImpl.h"

ParallelUnionFindImpl::ParallelUnionFindImpl(void)
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
