#include "ParallelUnionFind2DStripes.h"

ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info) :
ParallelUnionFindImpl(info),
wuf(new WeightedUnionFind(info.domainWidth * info.domainHeight))
{
    if(mDecompositionInfo.numOfProc <= 0)
    {
        std::cerr << "0 number of processors!" << std::endl;
        std::cerr << "Check whether the MPI has been initialized!" << std::endl;
    }
    else
    {
        std::cout << "Processor " << mDecompositionInfo.myRank << " of " << mDecompositionInfo.numOfProc << std::endl;
        std::cout << "PUF of type 2DStripes created, width " 
            << mDecompositionInfo.domainWidth << " height " << mDecompositionInfo.domainHeight << std::endl;
    }
}

ParallelUnionFind2DStripes::~ParallelUnionFind2DStripes(void)
{
}

void ParallelUnionFind2DStripes::runLocalUnionFind(void)
{
}

void ParallelUnionFind2DStripes::constructGlobalLabeling(void)
{
}

void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
}

void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
}
